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
 * \brief @copybrief Dumux::FluidSystems::H2ON2O2
 */
#ifndef DUMUX_H2O_N2_O2_FLUID_SYSTEM_HH
#define DUMUX_H2O_N2_O2_FLUID_SYSTEM_HH

#include <cassert>

#include <dumux/material/idealgas.hh>
#include <dumux/material/constants.hh>

#include <dumux/material/components/n2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/o2.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
#include <dumux/material/binarycoefficients/n2_o2.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>

#ifdef DUMUX_PROPERTIES_HH
#include <dumux/common/basicproperties.hh>
#include <dumux/material/fluidsystems/defaultcomponents.hh>
#endif

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase (water and air) fluid system
 *        with water, nitrogen and oxygen as components.
 *
 * This fluidsystem uses tabulated version of water of the IAPWS-formulation.
 *
 * Also remember to initialize tabulated components (FluidSystem::init()), while this
 * is not necessary for non-tabularized ones.
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace FluidSystem::.
 * An adapter class using FluidSystem<TypeTag> is also provided
 * at the end of this file.
 */
template <class Scalar, bool useComplexRelations = true>
class H2ON2O2
    : public BaseFluidSystem<Scalar, H2ON2O2<Scalar, useComplexRelations> >
{
    typedef H2ON2O2<Scalar, useComplexRelations> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    typedef Dumux::IdealGas<Scalar> IdealGas;
    typedef Dumux::Constants<Scalar> Constants;
    typedef Dumux::H2O<Scalar> IapwsH2O;
    typedef TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;
    typedef Dumux::N2<Scalar> SimpleN2;
    typedef Dumux::O2<Scalar> O2;

    //! The components for pure water
    typedef TabulatedH2O H2O;

    //! The components for pure nitrogen
    typedef SimpleN2 N2;

public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;
    //! Number of solid phases besides the solid matrix
    static constexpr int numSPhases = 3;   // Admin Bcmin Amin mineraux

    static constexpr int wPhaseIdx = 0; // index of the wetting phase
    static constexpr int nPhaseIdx = 1; // index of the non-wetting phase

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static const std::string name[] = {
            std::string("l"),
            std::string("g")
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != nPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if Henry's law and Raoult's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Raoult's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        if (phaseIdx == nPhaseIdx)
            return true;
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        if (phaseIdx == nPhaseIdx)
            // let the components decide
            return H2O::gasIsIdeal() && N2::gasIsIdeal() && O2::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 9;
    // There are no secondary components (still 2pNc model needs this parameter)
    static const int numSecComponents = 7;
    static const int numMajorComponents = 2;

    static constexpr int H2OIdx = 0;//first major component
    static constexpr int N2Idx = 1;//second major component

    static constexpr int HIdx   = 2;//
    static constexpr int AaqIdx = 3;//
    static constexpr int BaqIdx = 4;//
    static constexpr int CaqIdx = 5;//
    static constexpr int DIdx = 6;//


    static constexpr int AminIdx = 7;//
    static constexpr int ADminIdx = 8;//

    
    

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20C:
    static const int wCompIdx = wPhaseIdx; //=0 -> H2OIdx
    static const int nCompIdx = nPhaseIdx; //=1 -> N2Idx

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
            {
            case H2OIdx: return H2O::name();
            case N2Idx: return N2::name();
            case HIdx: return "H+";      
            case AaqIdx: return "Aaq";       
            case BaqIdx: return "Baq";   
            case CaqIdx: return "Caq";   
            case DIdx: return "D-";   
            case AminIdx: return "Amin";                              
            case ADminIdx: return "ADmin";                              
            }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx) // KG PAR MOL
    {
        static const Scalar M[] = {
            18.0153e-3, // H2O::molarMass(),
            N2::molarMass(),
            1.0078e-3, // H+
            100e-3,   // Aaq
            50e-3,     // Baq
            120e-3,    // Caq
            40e-3,     // D-
            
        };

        assert(0 <= compIdx && compIdx < numComponents);
        //std::cout<<"numcomponents = "<<numComponents<<std::endl;
        //std::cout<<"N2 molar mass"<< N2::molarMass()<<std::endl;

        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            H2O::criticalTemperature(),
            N2::criticalTemperature(),
            O2::criticalTemperature()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            H2O::criticalPressure(),
            N2::criticalPressure(),
            O2::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }

    /*!
     * \brief Molar volume of a component at the critical point \f$\mathrm{[m^3/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "H2ON2O2FluidSystem::criticalMolarVolume()");
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(),
            N2::acentricFactor(),
            O2::acentricFactor()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    }

    /*!
     * \brief Kelvin equation in \f$\mathrm{[Pa]}\f$
     *
     * Calculate the increase vapor pressure over the
     * curved surface of a drop with radius r
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     * \param radius The radius of the drop
     */
    template <class FluidState>
    static Scalar kelvinVaporPressure(const FluidState &fluidState,
                                      const int phaseIdx,
                                      const int compIdx,
                                      const Scalar radius)
    {
        assert(0 <= phaseIdx  && phaseIdx == wPhaseIdx);
        assert(0 <= compIdx  && compIdx == wCompIdx);

        Scalar T = fluidState.temperature(phaseIdx);

        Scalar vaporPressure = H2O::vaporPressure(T);
        Scalar exponent = molarMass(compIdx)/(density(fluidState, phaseIdx) * Constants::R * T);
        exponent *= (2 * surfaceTension(fluidState) / radius);
        using std::exp;
        Scalar kelvinVaporPressure = vaporPressure * exp(exponent);

        return kelvinVaporPressure;
    }

    /*!
     * \brief Calculate the surface tension between water and air in \f$\mathrm{[\frac{N}{m}]}\f$,
     * according to IAPWS Release on Surface Tension from September 1994.
     * The equation is valid between the triple Point (0.01C) and the critical temperature.
     *
     * \param fluidState An abitrary fluid state
     */
    template <class FluidState>
    static Scalar surfaceTension(const FluidState &fluidState)
    {
        const Scalar T = fluidState.temperature(); //K
        const Scalar B   = 0.2358 ; // [N/m]
        const Scalar T_c = H2O::criticalTemperature(); //K
        const Scalar mu  = 1.256;
        const Scalar b   = -0.625;
        //Equation to calculate surface Tension of Water According to IAPWS Release on Surface Tension from September 1994
        using std::pow;
        const Scalar surfaceTension = B*pow((1.-(T/T_c)),mu)*(1.+b*(1.-(T/T_c)));
        return surfaceTension; //surface Tension [N/m]
    }
    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (useComplexRelations)
            std::cout << "Using complex H2O-N2-O2 fluid system\n";
        else
            std::cout << "Using fast H2O-N2-O2 fluid system\n";

        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * If useComplexRelations == true, we apply Eq. (7)
     * in Class et al. (2002a) \cite A3:class:2002b <BR>
     * for the liquid density.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
                // assume pure water
            return 997.0 ; // change : given for the water phase       
        }

        // gas phase
        using std::max;
            // for the gas phase assume an ideal gas
            return
                1.13 ; // change : given as the density of the Gas
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * Compositional effects in the gas phase are accounted by the Wilke method.
     * See Reid et al. (1987)  \cite reid1987 <BR>
     * 4th edition, McGraw-Hill, 1987, 407-410
     * 5th edition, McGraw-Hill, 20001, p. 9.21/22
     * \note Compositional effects for a liquid mixture have to be implemented.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);


        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            // assume pure water for the liquid phase
            return 8.9e-4; // change water viscosity
        }

        // gas phase
        // assume pure nitrogen for the gas phase
        return 1.78e-5; // change gaz viscosity
    }



    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * For liquids with very low miscibility this boils down to the
     * inverse Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            switch(compIdx){
            case H2OIdx: return H2O::vaporPressure(T)/p;
            case N2Idx: return BinaryCoeff::H2O_N2::henry(T)/p;
            };
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\kappa\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mu_\kappa\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

#ifndef NDEBUG
        if (compIIdx == compJIdx ||
            phaseIdx > numPhases - 1 ||
            compJIdx > numComponents - 1)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
#endif

 

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            return 8e-11; //  // change Water diffusion
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
        // gas phase
        if (phaseIdx == nPhaseIdx) {
            return 3e-6; // change Gas diffusion
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }

        DUNE_THROW(Dune::InvalidStateException,
                  "Binary diffusion coefficient of components "
                  << compIIdx << " and " << compJIdx
                  << " in phase " << phaseIdx << " is undefined!\n");
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     *  \note This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        // gas phase
        else if (phaseIdx == nPhaseIdx)
        {
            // assume ideal mixture: which means
            // that the total specific enthalpy is the sum of the
            // "partial specific enthalpies" of the components.
            Scalar hH2O =
                fluidState.massFraction(nPhaseIdx, H2OIdx)
                * H2O::gasEnthalpy(T, p);
            Scalar hN2 =
                fluidState.massFraction(nPhaseIdx, N2Idx)
                * N2::gasEnthalpy(T,p);
            return hH2O + hN2;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
    /*!
     * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
     */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int componentIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Component enthalpies");
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of air and water as a first approximation.
     *
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar temperature  = fluidState.temperature(phaseIdx) ;
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else
        {
            Scalar lambdaPureN2 = N2::gasThermalConductivity(temperature, pressure);
            Scalar lambdaPureO2 = O2::gasThermalConductivity(temperature, pressure);
            if (useComplexRelations)
            {
                Scalar xN2 = fluidState.moleFraction(phaseIdx, N2Idx);
                Scalar xH2O = fluidState.moleFraction(phaseIdx, H2OIdx);
                Scalar lambdaN2 = xN2 * lambdaPureN2;
                Scalar partialPressure  = pressure * xH2O;
                Scalar lambdaH2O = xH2O * H2O::gasThermalConductivity(temperature, partialPressure);
                return lambdaN2 + lambdaH2O;
            }
            else
                return lambdaPureN2;
        }
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg*K]}\f$.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }

        Scalar c_pN2;
        Scalar c_pH2O;
        // let the water and nitrogen components do things their own way

        // assume an ideal gas for both components. See:
        //
            //http://en.wikipedia.org/wiki/Heat_capacity
        Scalar c_vN2molar = Constants::R*2.39;
        Scalar c_pN2molar = Constants::R + c_vN2molar;
        
        
        Scalar c_vH2Omolar = Constants::R*3.37; // <- correct??
        Scalar c_pH2Omolar = Constants::R + c_vH2Omolar;
        
        c_pN2 = c_pN2molar/molarMass(N2Idx);
        c_pH2O = c_pH2Omolar/molarMass(H2OIdx);
        
        
        // mangle all components together
        return
            c_pH2O*fluidState.massFraction(nPhaseIdx, H2OIdx)
            + c_pN2*fluidState.massFraction(nPhaseIdx, N2Idx);
    }

};

} // end namespace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
/*!
 * \brief A two-phase (water and air) fluid system
 *        with water, nitrogen and oxygen as components.
 *
 * This is an adapter to use H2ON2O2<TypeTag>, as is
 * done with most other classes in Dumux.
 */
template<class TypeTag>
class H2ON2O2FluidSystem
: public FluidSystems::H2ON2O2<typename GET_PROP_TYPE(TypeTag, Scalar),
                             GET_PROP_VALUE(TypeTag, EnableComplicatedFluidSystem)>
{};
#endif

} // end namespace

#endif
