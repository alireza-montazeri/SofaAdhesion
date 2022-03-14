/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: Alireza Montazeri                                                  *
*                                                                             *
* Contact information: alireza.montazeri9675@gmail.com                        *
******************************************************************************/

#define SOFA_ADHESIONPLUGIN_CPP

#include <AdhesionPlugin.inl>
#include <sofa/core/ObjectFactory.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa::component::forcefield
{

using namespace sofa::defaulttype;


int AdhesionPluginClass = core::RegisterObject("Adhesion Plugin")
        .add<AdhesionPlugin<Vec3Types> >(true);

template class SOFAADHESION_API AdhesionPlugin<Vec3Types>;

} // namespace sofa::component::forcefield
