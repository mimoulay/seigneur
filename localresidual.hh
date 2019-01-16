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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase n-component box model.
 */

#ifndef DUMUX_2PNC_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_2PNC_LOCAL_RESIDUAL_BASE_HH

#include "properties.hh"
#include "volumevariables.hh"
#include <dumux/nonlinear/newtoncontroller.hh>

#include <iostream>
#include <vector>

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase n-component fully implicit box model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the two-phase n-component flow.
 */
template<class TypeTag>
class TwoPNCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef TwoPNCLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        conti0EqIdx = Indices::conti0EqIdx,
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
        
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwsn = TwoPNCFormulation::pwsn,
        pnsw = TwoPNCFormulation::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::ctype CoordScalar;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPNCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the storage term of the current solution in a
     *        single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void evalPhaseStorage(const Element &element, int phaseIdx)
    {
        FVElementGeometry fvGeometry;
        fvGeometry.update(this->gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(this->problem_(), element, fvGeometry);
        ElementVolumeVariables volVars;
        volVars.update(this->problem_(), element, fvGeometry, false);

        this->residual_.resize(fvGeometry.numScv);
        this->residual_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeometry;
        this->bcTypesPtr_ = &bcTypes;
        this->prevVolVarsPtr_ = 0;
        this->curVolVarsPtr_ = &volVars;
        evalPhaseStorage_(phaseIdx);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage the mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const auto& elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const auto& volVars = elemVolVars[scvIdx];

        // Compute storage term of all fluid components in the fluid phases
        storage = 0;

        Scalar phi = volVars.porosity();

        Scalar Sl  = volVars.saturation(wPhaseIdx);
        Scalar Sg  = volVars.saturation(nPhaseIdx);
            
        Scalar rhoW= volVars.molarDensity(wPhaseIdx);
        Scalar rhoN= volVars.molarDensity(nPhaseIdx);

        //Major components :

        Scalar H2O = volVars.moleFraction(wPhaseIdx,  wCompIdx);
        Scalar N2l = volVars.moleFraction(wPhaseIdx,  nCompIdx);

        Scalar N2g = volVars.moleFraction(nPhaseIdx,  nCompIdx);

        // Primary components
        // mole fraction components of aqeuous phase
      
        Scalar H = volVars.moleFraction(wPhaseIdx, HIdx);
        Scalar Aaq = volVars.moleFraction(wPhaseIdx, AaqIdx);
        Scalar Baq = volVars.moleFraction(wPhaseIdx, BaqIdx);
        Scalar Caq = volVars.moleFraction(wPhaseIdx, CaqIdx);
        Scalar D = volVars.moleFraction(wPhaseIdx, DIdx);

        // Mineraux 
        Scalar Amin = volVars.moleFraction(wPhaseIdx, AminIdx);
        Scalar ADmin = volVars.moleFraction(wPhaseIdx, ADminIdx);


        // Reaction :  OH-  = H2O – H+ 
        Scalar aH= H/0.018;           // 0.18 water molality
        Scalar KOH=pow(10,-14);
        Scalar aOH=KOH/aH;        // aH2O = 1 
        Scalar OH=aOH*0.018;

        // Reaction : Dg = D-   + H+  -  H2O
        Scalar aD = D/0.018 ; 
        Scalar KDg = pow(10,5); // ou bien -5 ????
        Scalar fDg = KDg * aD * aH ;  // fugacity du gaz Dg
        Scalar Dg = fDg / volVars.pressure(nPhaseIdx) ; // Loi d'henry xDg* Pg= f*Dg
        
        // Majour components ; 
        storage[wCompIdx] += phi*rhoW*Sl*(H2O+ OH) - phi*rhoN*Dg
            + 2*phi*Amin ;

        storage[wCompIdx] += phi*rhoW*Sl*N2l + phi*rhoN*Sg*N2g;
            

        // RTM components
        
        storage[HIdx] += phi*rhoW*Sl*(H-OH) + phi*rhoN*Dg + phi* Amin ; // (1-phi)
        storage[AaqIdx] +=  phi*rhoW*Sl*Aaq + phi*(Amin+ ADmin) ; // (1-phi)
        storage[BaqIdx] +=  phi*rhoW*Sl*Baq ; // (1-phi)
        storage[CaqIdx] +=  phi*rhoW*Sl*Caq ;  // (1-phi)
        storage[DIdx] +=  phi*rhoW*Sl*D +  phi*rhoN*Sg*Dg + phi * ADmin; 
        

        Valgrind::CheckDefined(storage);
    }
    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fIdx The index of the sub-control-volume face
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, bool onBoundary = false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->curVolVars_(),
                        onBoundary);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
    }







    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
           
        Scalar qW  =  fluxVars.volumeFlux(wPhaseIdx);
        Scalar qN  =  fluxVars.volumeFlux(nPhaseIdx);
        
        const VolumeVariables &upW =
            this->curVolVars_(fluxVars.upstreamIdx(wPhaseIdx));
        const VolumeVariables &dnW =
            this->curVolVars_(fluxVars.downstreamIdx(wPhaseIdx));
        
        const VolumeVariables &upN =
            this->curVolVars_(fluxVars.upstreamIdx(nPhaseIdx));
        const VolumeVariables &dnN =
            this->curVolVars_(fluxVars.downstreamIdx(nPhaseIdx));
        
        if (massUpwindWeight_ > 0.0)
            {
                // upstream vertex
                Scalar rhoW =  upW.molarDensity(wPhaseIdx);
                Scalar rhoN =  upN.molarDensity(nPhaseIdx);
                
                Scalar H2O  =  upW.moleFraction(wPhaseIdx, wCompIdx);
                Scalar N2l  =  upW.moleFraction(wPhaseIdx, nCompIdx);
                Scalar N2g  =  upN.moleFraction(nPhaseIdx, nCompIdx);
            
                Scalar H    =  upW.moleFraction(wPhaseIdx, HIdx  );
                Scalar Aaq  =  upW.moleFraction(wPhaseIdx, AaqIdx );
                Scalar Baq  =  upW.moleFraction(wPhaseIdx, BaqIdx );
                Scalar Caq  =  upW.moleFraction(wPhaseIdx, CaqIdx );
                Scalar D    =  upW.moleFraction(wPhaseIdx, DIdx  );
                     
                // Reaction :  OH-  = H2O – H+ 
                Scalar aH   =  H/0.018;           // 0.18 water molality
                Scalar KOH  =  pow(10,-14);
                Scalar aOH  =  KOH/aH;        // aH2O = 1 
                Scalar OH   =  aOH*0.018;
                
                // Reaction : Dg = D-   + H+  -  H2O
                Scalar aD   = D/0.018 ; 
                Scalar KDg  = pow(10,5); // ou bien -5 ????
                Scalar fDg  = KDg * aD * aH ;  // fugacity du gaz Dg
                Scalar Dg   = fDg /  upN.pressure(nPhaseIdx) ; // Loi d'henry xDg* Pg= f*Dg
                
                // massUpwindWeight
                flux[wCompIdx] += (qW*rhoW* (H2O+OH) - qN*rhoN* Dg   ) * massUpwindWeight_ ;
                flux[nCompIdx] += (qW*rhoW* N2l      + qN*rhoN* N2g  ) * massUpwindWeight_ ;
                flux[HIdx]     += (qW*rhoW* (H -OH)  + qN*rhoN* Dg   ) * massUpwindWeight_ ;
                flux[AaqIdx]   += (qW*rhoW* Aaq                      ) * massUpwindWeight_ ;
                flux[BaqIdx]   += (qW*rhoW* Baq                      ) * massUpwindWeight_ ;
                flux[CaqIdx]   += (qW*rhoW* Caq                      ) * massUpwindWeight_ ;
                flux[DIdx]     += (qW*rhoW* D   + qN*rhoN* Dg        ) * massUpwindWeight_ ;
             

            }
        
        if (massUpwindWeight_ < 1.0)
            {
                // downstream vertex
                Scalar rhoW =  dnW.molarDensity(wPhaseIdx);
                Scalar rhoN =  dnW.molarDensity(nPhaseIdx);
                
                Scalar H2O  =  dnW.moleFraction(wPhaseIdx, wCompIdx);
                Scalar N2l  =  dnW.moleFraction(wPhaseIdx, nCompIdx);
                Scalar N2g  =  dnN.moleFraction(nPhaseIdx, nCompIdx);
            
                Scalar H    =  dnW.moleFraction(wPhaseIdx, HIdx  );
                Scalar Aaq  =  dnW.moleFraction(wPhaseIdx, AaqIdx );
                Scalar Baq  =  dnW.moleFraction(wPhaseIdx, BaqIdx );
                Scalar Caq  =  dnW.moleFraction(wPhaseIdx, CaqIdx );
                Scalar D    =  dnW.moleFraction(wPhaseIdx, DIdx  );
                     
                // Reaction :  OH-  = H2O – H+ 
                Scalar aH   =  H/0.018;           // 0.18 water molality
                Scalar KOH  =  pow(10,-14);
                Scalar aOH  =  KOH/aH;        // aH2O = 1 
                Scalar OH   =  aOH*0.018;
                
                // Reaction : Dg = D-   + H+  -  H2O
                Scalar aD   = D/0.018 ; 
                Scalar KDg  = pow(10,5); // ou bien -5 ????
                Scalar fDg  = KDg * aD * aH ;  // fugacity du gaz Dg
                Scalar Dg   = fDg /  dnN.pressure(nPhaseIdx) ; // Loi d'henry xDg* Pg= f*Dg
                
                // massUpwindWeight
                flux[wCompIdx] += (qW*rhoW* (H2O+OH) - qN*rhoN* Dg   ) * (1-massUpwindWeight_) ;
                flux[nCompIdx] += (qW*rhoW* N2l      + qN*rhoN* N2g  ) * (1-massUpwindWeight_) ;
                flux[HIdx]     += (qW*rhoW* (H -OH)  + qN*rhoN* Dg   ) * (1-massUpwindWeight_) ;
                flux[AaqIdx]   += (qW*rhoW* Aaq                      ) * (1-massUpwindWeight_) ;
                flux[BaqIdx]   += (qW*rhoW* Baq                      ) * (1-massUpwindWeight_) ;
                flux[CaqIdx]   += (qW*rhoW* Caq                      ) * (1-massUpwindWeight_) ;
                flux[DIdx]     += (qW*rhoW* D   + qN*rhoN* Dg        ) * (1-massUpwindWeight_) ;

            }

    }





    
    /*!
     * \brief Evaluates the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive flux of H2O, N2L ,H and Aaq Baq Caq and D component in liquid phase
        Scalar D_H2Ol= -(fluxVars.moleFractionGradW(wCompIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,wCompIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[wCompIdx] += D_H2Ol;
        
        Scalar D_N2l= -(fluxVars.moleFractionGradW(nCompIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,nCompIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[nCompIdx] += D_N2l;
            
        Scalar D_H= -(fluxVars.moleFractionGradW(HIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,HIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[HIdx] +=  D_H;
        
        Scalar D_Aaq= -(fluxVars.moleFractionGradW(AaqIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,AaqIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[AaqIdx] +=  D_Aaq;

        Scalar D_Baq= -(fluxVars.moleFractionGradW(BaqIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,BaqIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[BaqIdx] +=  D_Baq;
        
        Scalar D_Caq= -(fluxVars.moleFractionGradW(CaqIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,CaqIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[CaqIdx] +=  D_Caq;


        Scalar D_D= -(fluxVars.moleFractionGradW(DIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(wPhaseIdx,DIdx)
            * fluxVars.molarDensity(wPhaseIdx);
        flux[DIdx] +=  D_D;

        // diffusive flux of N2 in non-wetting phase        
        Scalar D_N2g = -(fluxVars.moleFractionGradN(nPhaseIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(nPhaseIdx,nCompIdx)
            * fluxVars.molarDensity(nPhaseIdx);
        flux[nCompIdx] += D_N2g;

        // contribution du Dg dans la diffusion des composant H20 H et D
        Scalar D_H2Og = -(fluxVars.moleFractionGradN(wCompIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(nPhaseIdx,wCompIdx)
            * fluxVars.molarDensity(nPhaseIdx);
        // Reaction Dg = -H2O + H + D 
        // le Dd contribue aux equations de H20 H ET D
        flux[wCompIdx] += D_H2Og;

        Scalar D_Hg = -(fluxVars.moleFractionGradN(HIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(nPhaseIdx,HIdx)
            * fluxVars.molarDensity(nPhaseIdx);
        flux[HIdx] += D_Hg;

        Scalar D_Dg = -(fluxVars.moleFractionGradN(DIdx)*fluxVars.face().normal)
            * fluxVars.porousDiffCoeff(nPhaseIdx,DIdx)
            * fluxVars.molarDensity(nPhaseIdx);
        flux[DIdx] += D_Dg;




    }


    void computeSource(PrimaryVariables &source, const int scvIdx)
    {
        
        source = 0;
        
        const ElementVolumeVariables &elemVolVars = this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];
        


        // Primary components
        // mole fraction components of aqeuous phase
      
        Scalar Aaq = volVars.moleFraction(wPhaseIdx, AaqIdx);
        Scalar D = volVars.moleFraction(wPhaseIdx, DIdx);

        // Mineraux 
        Scalar Amin = volVars.moleFraction(wPhaseIdx, AminIdx);
        Scalar ADmin = volVars.moleFraction(wPhaseIdx, ADminIdx);

        // constantes de reaction de precipitation
        Scalar KAmin = 10;
        Scalar KADmin = pow(10,12); 

        // les activités des primaires
        Scalar aAaq=Aaq/0.018;
        Scalar aD=D/0.018;
        
        // loi d'action de masse des minereaux 
        Scalar QAmin=KAmin*aAaq;
        Scalar QADmin=KADmin*aAaq*aD;

        source[AminIdx]  =std::min(Amin,1.-QAmin); 
        source[ADminIdx]  =std::min(ADmin,1.-QADmin); 
    
        
    }


protected:

    void evalPhaseStorage_(int phaseIdx)
    {
        // evaluate the storage terms of a single phase
        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            PrimaryVariables &result = this->residual_[scvIdx];
            const ElementVolumeVariables &elemVolVars = this->curVolVars_();
            const VolumeVariables &volVars = elemVolVars[scvIdx];

            // compute storage term of all fluid components within all phases
            result = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                result[conti0EqIdx + compIdx] += volVars.density(phaseIdx)
                                                    * volVars.saturation(phaseIdx)
                                                    * volVars.massFraction(phaseIdx, compIdx)
                                                    * volVars.porosity();
            }
            result *= this->fvGeometry_().subContVol[scvIdx].volume;
        }
    }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

public:
   Scalar massUpwindWeight_;
};

} // end namespace

#endif
