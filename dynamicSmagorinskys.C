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

    
\*---------------------------------------------------------------------------*/

#include "dynamicSmagorinsky.H"

#include "turbulentFluidThermoModels.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



makeLESModel(dynamicSmagorinsky);
                                            


// ************************************************************************* //
