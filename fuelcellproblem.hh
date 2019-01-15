// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
#ifndef DUMUX_FUELCELL_PROBLEM_HH
#define DUMUX_FUELCELL_PROBLEM_HH

#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>

//#include <dumux/porousmediumflow/2pnc/implicit/model.hh>
//#include <dumux/material/fluidsystems/h2on2o2.hh>

#include "model/model.hh"
#include "fluidsystem.hh"
#include <dumux/linear/amgbackend.hh>


#include "fuelcellspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class FuelCellProblem;

namespace Properties
{
NEW_TYPE_TAG(FuelCellProblem, INHERITS_FROM(TwoPNC, FuelCellSpatialParams));
NEW_TYPE_TAG(FuelCellBoxProblem, INHERITS_FROM(BoxModel, FuelCellProblem));
NEW_TYPE_TAG(FuelCellCCProblem, INHERITS_FROM(CCModel, FuelCellProblem));

// Set the grid type
// change from 2d to 1d
SET_TYPE_PROP(FuelCellProblem, Grid, Dune::YaspGrid<1>);
// Set the problem property
SET_TYPE_PROP(FuelCellProblem, Problem, FuelCellProblem<TypeTag>);
// Set the primary variable combination for the 2pnc model
SET_INT_PROP(FuelCellProblem, Formulation, TwoPNCFormulation::pnsw);

// Set fluid configuration
SET_PROP(FuelCellProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // change pour simplifier les calculs densité et autres dans fluidsystem
    //static const bool useComplexRelations = true;
    static const bool useComplexRelations = false;

 public:
    typedef FluidSystems::H2ON2O2<Scalar, useComplexRelations> type;
};

// Set the transport equation that is replaced by the total mass balance
SET_INT_PROP(FuelCellProblem, ReplaceCompEqIdx, 3);
}

    // change
// use the algebraic multigrid
//SET_TYPE_PROP(FuelCellProblem, LinearSolver, AMGBackend<TypeTag> );

/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class FuelCellProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        numComponents = FluidSystem::numComponents,
        numSecComponents = FluidSystem::numSecComponents,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        wCompIdx = FluidSystem::wCompIdx, //major component of the liquid phase
        nCompIdx = FluidSystem::nCompIdx, //major component of the gas phase
        HIdx=2,
        AaqIdx=3,
        BaqIdx=4,
        CaqIdx=5,
        DIdx=6,
        AminIdx=7,
        BCminIdx=8,
        ADminIdx=9,
    };
    enum
    {
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        switchIdx = Indices::switchIdx, //liquid saturation or mole fraction
        conti0EqIdx = Indices::conti0EqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    // Select the electrochemistry method
    typedef Dumux::ElectroChemistry<TypeTag, ElectroChemistryModel::Ochs> ElectroChemistry;
    typedef Constants<Scalar> Constant;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FuelCellProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        nTemperature_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NTemperature);
        nPressure_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NPressure);
        pressureLow_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow);
        pressureHigh_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh);
        temperatureLow_     = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow);
        temperatureHigh_    = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh);
        temperature_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, InitialTemperature);

        name_               = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        pO2Inlet_            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, ElectroChemistry, pO2Inlet);

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \name Problem parameters
     */

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

   
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < eps_ )
            values.setAllDirichlet();
        else 
            values.setAllNeumann();
    }

  
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {initial_(values, globalPos);}


    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    { priVars = 0; }

    
   
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a sub control volume.
     *
     * \param element The element of the sub control volume
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The sub control volume index
     */

    int initialPhasePresence(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             int scvIdx) const
    {
        return Indices::bothPhases;
    }

  

private:

    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
       Scalar x=globalPos[0];
        
        if(x<2.0+eps_)
            {
                values[wPhaseIdx] =0.6;// 1e5/(997*9.81);
                values[nPhaseIdx] =1e5;// 1e5/(997*9.81);

                values[HIdx]= std::log(0.018*pow(10,-5)); //H+
                values[AaqIdx]=  std::log(1e-20); //Aaq
                values[BaqIdx]= std::log(0.3*0.018); ; //Baq
                values[CaqIdx]=  std::log(1e-20);//Caq
                values[DIdx]=   std::log(0.08211959342457908);//D-

                values[AminIdx]=   0.;//Amin
                values[ADminIdx]=  0.; //ADmin 
            }
        else if (x>=2.0 && x<=52.0+ eps_)
            {
                values[wPhaseIdx] =0.4;// 1e5/(997*9.81);;
                values[nPhaseIdx] =1e5;// 1e5/(997*9.81);

                values[HIdx]= std::log(0.018*pow(10,-5)); //H+
                values[AaqIdx]=  std::log(0.1*0.018); //Aaq
                values[BaqIdx]= std::log(1e-20); ; //Baq
                values[CaqIdx]=  std::log(0.1);//Caq
                values[DIdx]=   std::log(1e-20);//D-

                values[AminIdx]=   1600;//Amin  mol/m^3
                values[ADminIdx]=  0.; //ADmin  
            }
      
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    bool inReactionLayer_(const GlobalPosition& globalPos) const
    { return globalPos[1] < 0.1*(this->bBoxMax()[1] - this->bBoxMin()[1]) + eps_; }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    int nTemperature_;
    int nPressure_;
    std::string name_ ;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar pO2Inlet_;
};

} //end namespace

#endif
