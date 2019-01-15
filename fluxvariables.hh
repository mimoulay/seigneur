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
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase two-component model fully implicit model.
 */
#ifndef DUMUX_2PNC_FLUX_VARIABLES_HH
#define DUMUX_2PNC_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase n-component fully implicit model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */

template <class TypeTag>
class TwoPNCFluxVariables : public GET_PROP_TYPE(TypeTag, BaseFluxVariables)
{
    friend typename GET_PROP_TYPE(TypeTag, BaseFluxVariables); // be friends with base class
    typedef typename GET_PROP_TYPE(TypeTag, BaseFluxVariables) BaseFluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld,
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
          };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
            wPhaseIdx = FluidSystem::wPhaseIdx,
            nPhaseIdx = FluidSystem::nPhaseIdx,
            wCompIdx  = FluidSystem::wCompIdx,
            nCompIdx  = FluidSystem::nCompIdx,

            // change 
            HIdx=2,
            AaqIdx=3,
            BaqIdx=4,
            CaqIdx=5,
            DIdx=6,
            AminIdx=7,
            ADminIdx=8,
            
            // For diffusion l'indice  du Dg qui en qlq sorte remplit
            // la place du H20g 
            DgIdx = wPhaseIdx,            
         };

public:
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        BaseFluxVariables::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);
        calculatePorousDiffCoeff_(problem, element, elemVolVars);

        // change 

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                density_[phaseIdx] = Scalar(0);
                molarDensity_[phaseIdx] = Scalar(0);
            }
        

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                moleFractionGradW_[compIdx] = Scalar(0);
                moleFractionGradN_[compIdx] = Scalar(0);
                
            }
        /*====================================================*/


    }

protected:

    void calculateIpDensities_(const Problem &problem,
                               const Element &element,
                               const ElementVolumeVariables &elemVolVars)
    {
        // calculate densities at the integration points of the face
        density_.fill(0.0);
        molarDensity_.fill(0.0);
        for (unsigned int idx = 0; idx < this->face().numFap; idx++) // loop over adjacent vertices
        {
            // index for the element volume variables
            int volVarsIdx = this->face().fapIndices[idx];

            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                density_[phaseIdx] += elemVolVars[volVarsIdx].density(phaseIdx)*this->face().shapeValue[idx];
                molarDensity_[phaseIdx] += elemVolVars[volVarsIdx].molarDensity(phaseIdx)*this->face().shapeValue[idx];
            }
        }
    }

    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        calculateIpDensities_(problem, element, elemVolVars);
        BaseFluxVariables::calculateGradients_(problem, element, elemVolVars);

        /*
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            moleFractionGrad_[phaseIdx].fill(GlobalPosition(0.0));
            }*/

         // loop over number of flux approximation points
        for (unsigned int idx = 0; idx < this->face().numFap; ++idx)
        {
            // FE gradient at vertex idx
            const GlobalPosition &feGrad = this->face().grad[idx];

            // index for the element volume variables
            auto volVarsIdx = this->face().fapIndices[idx];

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            /*
            for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if(compIdx != phaseIdx) //No grad is needed for this case
                    {
                        moleFractionGrad_[phaseIdx][compIdx].axpy(elemVolVars[volVarsIdx].moleFraction(phaseIdx, compIdx), feGrad);
                    }
                }
           } 
            */

            // aqueous phase
            Scalar H2O =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, wCompIdx);
            Scalar N2l =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, nCompIdx);
            Scalar H    =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, HIdx);
            Scalar Aaq   =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, AaqIdx);
            Scalar Baq   =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, BaqIdx);
            Scalar Caq   =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, CaqIdx);
            Scalar D    =  elemVolVars[volVarsIdx].moleFraction(wPhaseIdx, DIdx);

            
            // gaz phase 
            Scalar N2g =  elemVolVars[volVarsIdx].moleFraction(nPhaseIdx,  nCompIdx);

            
            // Reaction :  OH-  = H2O – H+ 
            Scalar aH= H/0.018;           // 0.18 water molality
            Scalar KOH=pow(10,-14);
            Scalar aOH=KOH/aH;        // aH2O = 1 
            Scalar OH=aOH*0.018;

            // Reaction : Dg = D-   + H+  -  H2O
            Scalar aD = D/0.018 ; 
            Scalar KDg = pow(10,5); // ou bien -5 ????
            Scalar fDg = KDg * aD * aH ;  // fugacity du gaz Dg
            Scalar Dg = fDg /  elemVolVars[volVarsIdx].pressure(nPhaseIdx) ; // Loi d'henry xDg* Pg= f*Dg


            // Total concentration pour la diffusion
            Scalar TwH2O =  H2O +  OH  ; // aqueous part 
            Scalar TgH2O = -Dg ;              
            
            Scalar TwN2 = N2l ;
            Scalar TgN2 = N2g ;

            Scalar TwH =  H -  OH; 
            Scalar TgH = Dg ;

            Scalar TwAaq = Aaq ;  
            Scalar TwBaq = Baq ;  
            Scalar TwCaq = Caq ;  


            Scalar TwD = D; 
            Scalar TgD = Dg; 


             // gradient Total concentration pour la diffusion
            GlobalPosition grad_TwH2O = feGrad;
            GlobalPosition grad_TgH2O = feGrad;              

            GlobalPosition grad_TwN2 = feGrad ;
            GlobalPosition grad_TgN2 = feGrad ;

            GlobalPosition grad_TwH  = feGrad;
            GlobalPosition grad_TgH = feGrad ;

            GlobalPosition grad_TwAaq  = feGrad;
            GlobalPosition grad_TwBaq  = feGrad;
            GlobalPosition grad_TwCaq  = feGrad;
            
            GlobalPosition grad_TwD  = feGrad;
            GlobalPosition grad_TgD = feGrad ;




            grad_TwH2O *= TwH2O;
            grad_TgH2O *= TgH2O;              
            
            grad_TwN2 *= TwN2 ;
            grad_TgN2 *= TgN2 ;
            
            grad_TwH  *= TwH ;
            grad_TgH *=  TgH;
            
            grad_TwAaq  *= TwAaq;
            grad_TwBaq  *= TwBaq;
            grad_TwCaq  *= TwCaq;
            
            grad_TwD  *= TwD;
            grad_TgD *= TgD;
            
            // water gradient
            moleFractionGradW_[wPhaseIdx]+= grad_TwH2O;
            moleFractionGradW_[nPhaseIdx]+= grad_TwN2;
            moleFractionGradW_[HIdx]+= grad_TwN2;
            moleFractionGradW_[AaqIdx]+= grad_TwAaq;
            moleFractionGradW_[BaqIdx]+= grad_TwBaq;
            moleFractionGradW_[CaqIdx]+= grad_TwCaq;
            moleFractionGradW_[DIdx]+= grad_TwD;

            
            // gaz gradient 
            moleFractionGradN_[wPhaseIdx]+= grad_TgH2O;
            moleFractionGradN_[nPhaseIdx]+= grad_TgN2;
            moleFractionGradN_[HIdx]+= grad_TgN2;
            moleFractionGradN_[DIdx]+= grad_TgD;



        }
    }

    void calculatePorousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[this->face().i];
        const VolumeVariables &volVarsJ = elemVolVars[this->face().j];

        // the effective diffusion coefficients at vertex i and j
        Scalar diffCoeffI;
        Scalar diffCoeffJ;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            /* If there is no phase saturation on either side of the face
                * no diffusion takes place */
            if (volVarsI.saturation(phaseIdx) <= 0 || volVarsJ.saturation(phaseIdx) <= 0)
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    porousDiffCoeff_[phaseIdx][compIdx] = 0.0;
                }
            }
            else
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if (phaseIdx == compIdx)
                    {
                        porousDiffCoeff_[phaseIdx][compIdx] = 0.0;
                    }
                    else
                    {
                        diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                                     volVarsI.saturation(phaseIdx),
                                                                                     volVarsI.diffCoeff(phaseIdx, compIdx));

                        diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                                     volVarsJ.saturation(phaseIdx),
                                                                                     volVarsJ.diffCoeff(phaseIdx, compIdx));

                        porousDiffCoeff_[phaseIdx][compIdx] = harmonicMean(diffCoeffI, diffCoeffJ);
                    }
                }
            }
        }
    }

public:
    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     *
     *   \param phaseIdx The phase index
     *   \param compIdx The component index
     */
    Scalar porousDiffCoeff(int phaseIdx, int compIdx) const
    { return porousDiffCoeff_[phaseIdx][compIdx];}

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*!
     * \brief The mole fraction gradient of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    // change 
    const GlobalPosition &moleFractionGradW( int compIdx) const
    { return moleFractionGradW_[compIdx]; }

    const GlobalPosition &moleFractionGradN(int compIdx) const
    { return moleFractionGradN_[compIdx]; }

protected:
    // change 
    // mole fraction gradient
    //std::array<std::array<GlobalPosition, numComponents>, numPhases> moleFractionGradW_;
    //std::array<std::array<GlobalPosition, numComponents>, numPhases> moleFractionGradN_;

    GlobalPosition moleFractionGradW_[numComponents];
    GlobalPosition moleFractionGradN_[numComponents];

    // density of each face at the integration point
    std::array<Scalar, numPhases> density_, molarDensity_;

    // the diffusion coefficient for the porous medium
    Dune::FieldMatrix<Scalar, numPhases, numComponents> porousDiffCoeff_;
};

} // end namespace Dumux

#endif
