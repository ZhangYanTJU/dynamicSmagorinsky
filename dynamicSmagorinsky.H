/*---------------------------------------------------------------------------*\
dynamicSmagorinsky - Implementation of the dynamic Smagorinsky SGS model
		     as proposed by Lilly (1992) for OpenFOAM

Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2014 Alberto Passalacqua 

License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::compressible::LESModels::dynLocalAverageSmagorinsky

Description
    The dynamic Smagorinsky model for compressible flows.
     From Subgrid-Scale Models for Compressible LES,Martin, Piomelli, & Candler
     Theoretical Comput. Fluid Dynamics 2000 13:361-376


SourceFiles
    dynamicSmagorinsky.C
    
Author
    Ayodeji Ojofeitimi, ayodejio@umd.edu

Notes
    Implementation of the dynamic Smagorinsky model with coefficients cD and
    cI computed as local average of their face values to avoid numerical 
    instabilities. 

    Negative values of the effective viscosity are removed by clipping it to
    zero (muSgs is clipped to -mu)

    The code is known to work with OpenFOAM 6

\*---------------------------------------------------------------------------*/

#ifndef dynamicSmagorinsky_H
#define dynamicSmagorinsky_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynamicSmagorinsky Declaration
\*---------------------------------------------------------------------------*/
template<class BasicTurbulenceModel>
class dynamicSmagorinsky
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private data

        volScalarField k_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);
	
        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField cD(const volSymmTensorField& D) const;
        volScalarField cI(const volSymmTensorField& D) const;
        
        //- Calculate Ce by filtering the velocity field U
        volScalarField Ce
        (
            const volSymmTensorField& D,
            const volScalarField& KK
        ) const;        

        // Disallow default bitwise copy construct and assignment
        dynamicSmagorinsky(const dynamicSmagorinsky&);
        dynamicSmagorinsky& operator=(const dynamicSmagorinsky&);


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;
    
    
    //- Runtime type information
    TypeName("dynamicSmagorinsky");

    // Constructors

        //- Construct from components
        dynamicSmagorinsky
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

    //- Destructor
    virtual ~dynamicSmagorinsky()
    {}


    // Member Functions

    
        //- Read model coefficients if they have changed
        virtual bool read();
        
        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }
        
        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }
        
        //- Update the SGS eddy viscosity
        virtual void correctNut();
        
        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
	
        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const
        {
            Info<<"haha, Sct is controled by me!!!"<<endl;
            return this->mut()/Sct_ + this->mu();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicSmagorinsky.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
