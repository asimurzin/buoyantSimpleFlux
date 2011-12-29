#!/usr/bin/env python

#---------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey Simurzin
##


#---------------------------------------------------------------------------
from Foam import man, ref


#---------------------------------------------------------------------------
def createFields( runTime, mesh, g ):
    
    ref.ext_Info() << "Reading thermophysical properties\n" << ref.nl
    pThermo = man.basicPsiThermo.New( mesh );

    rho = man.volScalarField( man.IOobject( ref.word( "rho" ), 
                                            ref.fileName( runTime.timeName() ),
                                            mesh, 
                                            ref.IOobject.NO_READ, 
                                            ref.IOobject.NO_WRITE ),
                              man( pThermo.rho(), man.Deps( pThermo ) ) );
    
    p = man( pThermo.p(), man.Deps( pThermo ) )
    h = man( pThermo.h(), man.Deps( pThermo ) )
    psi = man( pThermo.psi(), man.Deps( pThermo ) )
    
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                              mesh )
    
    phi = man.compressibleCreatePhi( runTime, mesh, rho, U );

    ref.ext_Info() << "Creating turbulence model\n" << ref.nl
    turbulence = man.compressible.RASModel.New( rho, U,  phi, pThermo );

    ref.ext_Info()<< "Calculating field g.h\n" << ref.nl
    
    gh = man.volScalarField( ref.word( "gh" ), man( g & mesh.C(), man.Deps( mesh ) ) )
    ghf = man.surfaceScalarField( ref.word( "ghf" ), man( g & mesh.Cf(), man.Deps( mesh ) ) )
    
    ref.ext_Info() << "Reading field p_rgh\n" << ref.nl
    p_rgh = man.volScalarField( man.IOobject( ref.word( "p_rgh" ),
                                              ref.fileName( runTime.timeName() ),
                                              mesh,
                                              ref.IOobject.MUST_READ,
                                              ref.IOobject.AUTO_WRITE ),
                                mesh );
    # Force p_rgh to be consistent with p
    p_rgh << p - rho * gh
    
    pRefCell = 0
    pRefValue = 0.0
    pRefCell, pRefValue = ref.setRefCell( p, p_rgh, mesh.solutionDict().subDict( ref.word( "SIMPLE" ) ), pRefCell, pRefValue );

    initialMass = ref.fvc.domainIntegrate( rho )
    totalVolume = mesh.V().ext_sum()
    
    return pThermo, rho, p, h, psi, U, phi, turbulence, gh, ghf, p_rgh, pRefCell, pRefValue, initialMass, totalVolume


#--------------------------------------------------------------------------------------
def fun_Ueqn( simple, mesh, rho, U, phi, turbulence, ghf, p_rgh ):
    UEqn = man.fvm.div( phi, U ) + man( turbulence.divDevRhoReff( U ), man.Deps( turbulence, U) )
    UEqn.relax()
 
    if simple.momentumPredictor():
     ref.solve( UEqn == man.fvc.reconstruct( ( ( - ghf * man.fvc.snGrad( rho ) - man.fvc.snGrad( p_rgh ) )
                                               * man.surfaceScalarField( mesh.magSf(), man.Deps( mesh ) ) ) ) )

    return UEqn

#--------------------------------------------------------------------------------------
def fun_hEqn( thermo, rho, p, h, phi, U, turbulence ):
    hEqn = ref.fvm.div( phi, h ) - ref.fvm.Sp( ref.fvc.div( phi ), h ) - ref.fvm.laplacian(turbulence.alphaEff(), h ) \
            == - ref.fvc.div( phi, 0.5 * U.magSqr(), ref.word( "div(phi,K)" ) )

    hEqn.relax()
    hEqn.solve()
    
    thermo.correct()
    pass
    
    
#--------------------------------------------------------------------------------------
def fun_pEqn( mesh, runTime, simple, thermo, rho, p, h, psi, U, phi, turbulence, \
                      gh, ghf, p_rgh, UEqn, pRefCell, pRefValue, cumulativeContErr, initialMass):
      
      rho << thermo.rho()
      rho.relax()

      rAU = 1.0 / UEqn.A()
      rhorAUf = ref.surfaceScalarField( ref.word( "(rho*(1|A(U)))" ), ref.fvc.interpolate( rho() * rAU ) )
      
      U << rAU * UEqn.H()
      
      phi << ref.fvc.interpolate( rho ) * ( ref.fvc.interpolate( U() ) & mesh.Sf() )
      
      closedVolume = ref.adjustPhi( phi, U, p_rgh )
      
      buoyancyPhi = rhorAUf * ghf * ref.fvc.snGrad( rho ) * mesh.magSf()
      phi -= buoyancyPhi

      while simple.correctNonOrthogonal():

          p_rghEqn = ( ref.fvm.laplacian( rhorAUf, p_rgh ) == ref.fvc.div( phi ) )
          p_rghEqn.setReference( pRefCell, ref.getRefCellValue( p_rgh, pRefCell ) )
          p_rghEqn.solve()

          if simple.finalNonOrthogonalIter():
              # Calculate the conservative fluxes
              phi -= p_rghEqn.flux()
              
              # Explicitly relax pressure for momentum corrector
              p_rgh.relax()
              
              # Correct the momentum source with the pressure gradient flux
              # calculated from the relaxed pressure
              U -= rAU * ref.fvc.reconstruct( ( buoyancyPhi + p_rghEqn.flux() ) / rhorAUf )
              U.correctBoundaryConditions()
              pass
          pass
  
      cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )
      
      p << p_rgh + rho * gh

      # For closed-volume cases adjust the pressure level
      # to obey overall mass continuity
      if closedVolume:
          p += ( initialMass - ref.fvc.domainIntegrate( psi * p ) ) / ref.fvc.domainIntegrate( psi ) 
          p_rgh << p - rho * gh
          pass
      
      rho << thermo.rho()
      rho.relax()
      
      ref.ext_Info() << "rho max/min : " <<  rho.ext_max().value() << " " << rho.ext_min().value() << ref.nl
      
      return cumulativeContErr

    
                      
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    g = ref.readGravitationalAcceleration( runTime, mesh )

    pThermo, rho, p, h, psi, U, phi, turbulence, gh, ghf, p_rgh, \
         pRefCell, pRefValue, initialMass, totalVolume = createFields( runTime, mesh, g )

    cumulativeContErr = ref.initContinuityErrs()
    simple = man.simpleControl( mesh )

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl

    while simple.loop():
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl

        # Pressure-velocity SIMPLE corrector
        UEqn = fun_Ueqn( simple, mesh, rho, U, phi, turbulence, ghf, p_rgh )
        fun_hEqn( pThermo, rho, p, h, phi, U, turbulence )
        fun_pEqn( mesh, runTime, simple, pThermo, rho, p, h, psi, U, phi, turbulence, 
                  gh, ghf, p_rgh, UEqn, pRefCell, pRefValue, cumulativeContErr, initialMass)

        turbulence.correct()

        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
                       << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
                       << ref.nl << ref.nl
        pass

    ref.ext_Info() << "End\n" << ref.nl

    import os
    return os.EX_OK

    
#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020100" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
    re.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.0 or higher \n "     
    


    
#--------------------------------------------------------------------------------------

