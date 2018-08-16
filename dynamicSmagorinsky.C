/*---------------------------------------------------------------------------*\
dynamicSmagorinsky - Implementation of the dynamic Smagorinsky
		     SGS model.
    
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

\*---------------------------------------------------------------------------*/

#include "dynamicSmagorinsky.H"
#include "bound.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::updateSubGridScaleFields
(
    const volSymmTensorField& Sij
)
{
    // The SGS viscosity is bounded so that nuEff cannot become negative.
    // Values are limited here, and not in nuEff, for consistency in stored
    // data and in submodels using nuSgs().
    // No warning message is printed when this limitation is applied.
    this->nut_ = max(cD(Sij)*sqr(this->delta())*sqrt(magSqr(Sij)), -this->nu());
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();    
}

template<class BasicTurbulenceModel>
volScalarField dynamicSmagorinsky<BasicTurbulenceModel>::cD
(
    const volSymmTensorField& Sij
) const
{
    const rhoField& rho = this->rho_;
    volSymmTensorField Dfilter = filter_(rho*Sij)/filter_(rho);
    volScalarField magSijfilter = sqrt(2.0*(Dfilter && Dfilter));
    volScalarField magSij = sqrt(2*(Sij && Sij));

//  volSymmTensorField Lij = filter_(rho*sqr(U())) - sqr(filter_(rho*U()))/filter_(rho);
    volSymmTensorField Lij = filter_(sqr(rho*this->U_)/rho) - sqr(filter_(rho*this->U_))/filter_(rho);
    volSymmTensorField Bij = -2*sqr(2*this->delta())*filter_(rho)*magSijfilter*dev(Dfilter);
    volSymmTensorField Aij = -2*sqr(this->delta())*rho*magSij*dev(Sij);
    volSymmTensorField Mij = Bij - filter_(Aij);
    volScalarField LijMij = dev(Lij) && Mij;    //- Corrected Eq.17 in Mart et al (2000)
    volScalarField MklMkl = Mij && Mij;

     dimensionedScalar MM_Min 
    (
        "MM_Min",
        dimensionSet(2, -2, -4, 0, 0, 0, 0),
        scalar(SMALL)
     );
     
    dimensionedScalar LM_Min 
    (
        "LM_Min",
        dimensionSet(2, -2, -4, 0, 0, 0, 0),
        scalar(0.0)
     );    
    
    // Performing local average on cell faces and return
    return max(fvc::average(LijMij),LM_Min)/max(fvc::average(MklMkl),MM_Min);     
}

template<class BasicTurbulenceModel>
volScalarField dynamicSmagorinsky<BasicTurbulenceModel>::cI
(
    const volSymmTensorField& ij
) const
{
    tmp<volScalarField> KK = 
	0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)));

    const volScalarField mm
    (
        sqr(this->delta())*(4*sqr(mag(filter_(ij))) - filter_(sqr(mag(ij))))
    );

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    return fvc::average(KK*mm)/mmmm;
}

template<class BasicTurbulenceModel>
volScalarField dynamicSmagorinsky<BasicTurbulenceModel>::Ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volScalarField Ce
    (
        filter_(this->nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /filter_(pow(KK, 1.5)/(2.0*this->delta()))
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ce) + Ce);
    return tfld();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicSmagorinsky<BasicTurbulenceModel>::dynamicSmagorinsky
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    //updateSubGridScaleFields(symm(fvc::grad(U)));

    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correctNut()
{
}

template<class BasicTurbulenceModel>
void dynamicSmagorinsky<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    updateSubGridScaleFields(symm(fvc::grad(this->U_)));

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    const volSymmTensorField D(dev(symm(fvc::grad(this->U_))));
    const volScalarField G(this->GName(), 2.0*nut*(fvc::grad(this->U_) && D));
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), this->U_)));
    volScalarField KK(0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_))));
    KK.max(dimensionedScalar("small", KK.dimensions(), small));
    
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Ce(D, KK)*alpha*rho*sqrt(k_)/this->delta(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);
}

template<class BasicTurbulenceModel>
bool dynamicSmagorinsky<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
