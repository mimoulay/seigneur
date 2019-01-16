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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, n-component model.
 */
#ifndef DUMUX_2PNC_VOLUME_VARIABLES_HH
#define DUMUX_2PNC_VOLUME_VARIABLES_HH

#include <iostream>
#include <vector>

#include <dumux/implicit/model.hh>
// #include <dumux/material/fluidstates/compositional.hh>

#include <dumux/common/math.hh>

#include "compositional.hh"
#include "properties.hh"
#include "indices.hh"
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/miscible2pnccomposition.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class TypeTag>
class TwoPNCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        dim = GridView::dimension,
        dimWorld=GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

        // formulations
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPNCFormulation::pwsn,
        pnsw = TwoPNCFormulation::pnsw,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // component indices
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        // phase presence enums
        nPhaseOnly = Indices::nPhaseOnly,
        wPhaseOnly = Indices::wPhaseOnly,
        bothPhases = Indices::bothPhases,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        HIdx   = 2,//
        AaqIdx = 3,//
        BaqIdx = 4,//
        CaqIdx = 5,//
        DIdx = 6,//
        AminIdx = 7,//
        ADminIdx = 8,//

    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename Grid::ctype CoordScalar;
    typedef Dumux::Miscible2pNCComposition<Scalar, FluidSystem> Miscible2pNCComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;


    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
public:

    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     * \param priVars The primary Variables
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_, isOldSol);


        // change
        // Porosité initiale avant le changement 
        initialporosity_ = problem.spatialParams().porosity(element, fvGeometry, scvIdx);
        // volume flux initial 
        initial_VF = problem.spatialParams().initialVF(element, fvGeometry, scvIdx);

        
        Scalar Conc_Amin=fluidState_.moleFraction(wPhaseIdx,AminIdx);
        Scalar Conc_ADmin=fluidState_.moleFraction(wPhaseIdx,ADminIdx);
        Scalar VF=Conc_Amin*0.068e-3 + Conc_ADmin*0.218e-3;

        // Dispersivity
        dispersivity_ = problem.spatialParams().dispersivity(element, fvGeometry, scvIdx);

        // porosity
        porosity_ = initialporosity_-(VF-initial_VF);
        
        Valgrind::CheckDefined(porosity_);

        // Permeability factor
        permeabilityFactor_  =(std::pow(porosity_,3)*std::pow(1-initialporosity_,2))/(std::pow(1-porosity_,2)*std::pow(initialporosity_,3));


        /////////////
        // calculate the remaining quantities
        /////////////
        const MaterialLawParams &materialParams = problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting-phase saturation as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));

            mobility_[phaseIdx] = kr / fluidState_.viscosity(phaseIdx);
            Valgrind::CheckDefined(mobility_[phaseIdx]);
            int compIIdx = phaseIdx;
            for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                // binary diffusion coefficients
                diffCoeff_[phaseIdx][compJIdx] = 0.0;
                if(compIIdx!= compJIdx)
                {
                    diffCoeff_[phaseIdx][compJIdx] =
                        FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                paramCache,
                                                                phaseIdx,
                                                                compIIdx,
                                                                compJIdx);
                }
                Valgrind::CheckDefined(diffCoeff_[phaseIdx][compJIdx]);
            }
        }

      

        asImp_().updateEnergy_(priVars, problem,element, fvGeometry, scvIdx, isOldSol);
    }

   /*!
    * \copydoc ImplicitModel::completeFluidState
    * \param isOldSol Specifies whether this is the previous solution or the current one
    * \param priVars The primary Variables
    */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   int scvIdx,
                                   FluidState& fluidState,
                                   bool isOldSol = false)

    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        fluidState.setTemperature(t);

        int dofIdxGlobal = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        int phasePresence = problem.model().phasePresence(dofIdxGlobal, isOldSol);

        /////////////
        // set the saturations
        /////////////

        Scalar Sg;
        if (phasePresence == nPhaseOnly)
            Sg = 1.0;
        else if (phasePresence == wPhaseOnly)
        {
            Sg = 0.0;
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == pwsn)
                Sg = priVars[switchIdx];
            else if (formulation == pnsw)
                Sg = 1.0 - priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        }
        fluidState.setSaturation(nPhaseIdx, Sg);
        fluidState.setSaturation(wPhaseIdx, 1.0 - Sg);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams &materialParams
            = problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        Scalar pc = MaterialLaw::pc(materialParams, 1 - Sg);

        // extract the pressures
        if (formulation == pwsn) {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            if (priVars[pressureIdx] + pc < 0.0)
                 DUNE_THROW(NumericalProblem,"Capillary pressure is too low");
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw) {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            // Here we check for (p_g - pc) in order to ensure that (p_l > 0)
            if (priVars[pressureIdx] - pc < 0.0)
            {
                std::cout<< "p_g: "<< priVars[pressureIdx]<<" Cap_press: "<< pc << std::endl;
                DUNE_THROW(NumericalProblem,"Capillary pressure is too high");
            }
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;
        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases)
        {

            /* change : Pour le fluidestate on a 4 etapes : 
             */ 

            // etape 1 : calcul molefraction pour les composants 
            // secondaires dans la phase  liquid
            fluidState.setMoleFraction(wPhaseIdx, HIdx,   std::exp( priVars[HIdx] ) );
            fluidState.setMoleFraction(wPhaseIdx, AaqIdx,   std::exp( priVars[AaqIdx] ) );
            fluidState.setMoleFraction(wPhaseIdx, BaqIdx,   std::exp( priVars[BaqIdx] ) );
            fluidState.setMoleFraction(wPhaseIdx, CaqIdx,   std::exp( priVars[CaqIdx] ) );
            fluidState.setMoleFraction(wPhaseIdx, DIdx,   std::exp( priVars[DIdx] ) );

            // etape 2 : Calcul de xH2O liquid : 1- sum mole frac liquid
            Scalar sumMolFracliquid_ =std::exp( priVars[HIdx] )
                + std::exp( priVars[AaqIdx] )
                + std::exp( priVars[BaqIdx] ) 
                + std::exp( priVars[CaqIdx] )
                + std::exp( priVars[DIdx] ) ;
            fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1- sumMolFracliquid_);

            // etape 3 : Calcul du xN2g = 1-Dg 
            // On calcule d'abord Dg Reaction : Dg = D-   + H+  -  H2O
            Scalar aD = std::exp( priVars[DIdx] )/0.018 ;
            Scalar aH = std::exp( priVars[HIdx] )/0.018 ;
            Scalar KDg = pow(10,5); // ou bien -5 ????
            Scalar fDg = KDg * aD * aH ;  // fugacity du gaz Dg
            Scalar Dg = fDg / priVars[nPhaseIdx] ; // Loi d'henry xDg* Pg= f*Dg
            // N2g = 1-Dg
            fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1- Dg);

            // etape 4 : Calcul du xN2l par l'équilibre des phases. 
            // aN2l = fugacity_N2g equilibre des phases :
            // Pour un gaz parfait fugacity = pression partielle
            Scalar xN2g = 1-Dg ; 
            // aN2l = xN2g*Pg
            Scalar aN2l = xN2g *  priVars[nPhaseIdx] ; 
            Scalar N2l = aN2l * 0.018 ; 
            fluidState.setMoleFraction(wPhaseIdx, nCompIdx, N2l);

            
        }
        else if (phasePresence == nPhaseOnly)
        {
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            moleFrac[wCompIdx] =  priVars[switchIdx];
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
            }

            Scalar sumMoleFracOtherComponents = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                sumMoleFracOtherComponents+=moleFrac[compIdx];
            }

            sumMoleFracOtherComponents += moleFrac[wCompIdx];
            moleFrac[nCompIdx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(nPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             nPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);

        }
        else if (phasePresence == wPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            // extract _mass_ fractions in the nonwetting phase
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
            }
            moleFrac[nCompIdx] = priVars[switchIdx];
            Scalar sumMoleFracNotWater = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                sumMoleFracNotWater+=moleFrac[compIdx];
            }
            sumMoleFracNotWater += moleFrac[nCompIdx];
            moleFrac[wCompIdx] = 1 -sumMoleFracNotWater;


            // convert mass to mole fractions and set the fluid state
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
            fluidState.setViscosity(phaseIdx, mu);

            Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(FluidSystem::nPhaseIdx) - fluidState_.pressure(FluidSystem::wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffCoeff(int phaseIdx, int compIdx) const
    { return diffCoeff_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the molarity of a component in the phase \f$ moles/m^3 \f$
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     */
     Scalar molarity(int phaseIdx, int compIdx) const
    { return this->fluidState_.molarity(phaseIdx, compIdx);}

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     { return this->fluidState_.massFraction(phaseIdx, compIdx); }

    // change 
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }


     /*!
      * \brief Returns the mole fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     { return this->fluidState_.moleFraction(phaseIdx, compIdx); }


    Scalar permeabilityFactor() const
    { return permeabilityFactor_; }
    
    
protected:
    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return 0;
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                        const Problem &problem,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx,
                        bool isOldSol)
    { };

    Scalar porosity_; //!< Effective porosity within the control volume
    Scalar mobility_[numPhases]; //!< Effective mobility within the control volume
    Scalar density_;
    FluidState fluidState_;
    Scalar theta_;
    Scalar molWtPhase_[numPhases];
    Dune::FieldMatrix<Scalar, numPhases, numComponents> diffCoeff_;


    Scalar permeabilityFactor_;
    Scalar initialporosity_;
    Scalar initialpermeability_;
    Scalar dispersivity_ ;
    
    Scalar initial_VF;
    Scalar VF;


private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
