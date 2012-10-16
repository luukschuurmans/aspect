/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#include <aspect/postprocess/visualization/specific_heat.h>
#include <aspect/simulator.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SpecificHeat<dim>::
      SpecificHeat ()
        :
        DataPostprocessorScalar<dim> ("specific_heat",
                                      update_values | update_q_points)
      {}



      template <int dim>
      void
      SpecificHeat<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == dim+2+this->n_compositional_fields(), ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double pressure    = uh[q][dim];
            const double temperature = uh[q][dim+1];
            std::vector<double> composition(this->n_compositional_fields());
            for (unsigned int i=0;i<this->n_compositional_fields();++i)
              composition[i] = uh[q][dim+2+i];

            computed_quantities[q](0) = this->get_material_model().specific_heat(temperature,
                                                                                 pressure,
                                                                                 composition,
                                                                                 evaluation_points[q]);
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SpecificHeat,
                                                  "specific heat",
                                                  "A visualization output object that generates output "
                                                  "for the specific heat $C_p$.")
    }
  }
}
