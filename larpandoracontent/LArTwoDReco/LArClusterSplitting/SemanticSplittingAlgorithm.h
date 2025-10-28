/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/SemanticSplittingAlgorithm.h
 *
 *  @brief  Header file for the semantic splitting algorithm class.
 */

#ifndef LAR_SEMANTIC_SPLITTING_ALGORITHM_H
#define LAR_SEMANTIC_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

namespace lar_content
{

class SemanticSplittingAlgorithm : public TwoDSlidingFitSplittingAlgorithm
{
public:
    SemanticSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;
    
    /**
     *  @brief  Use sliding linear fit to identify the best split position
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the best split position based on semantic labels
     *
     *  @return pandora::StatusCode
     */
    pandora::StatusCode FindBestSplitPosition(
        const TwoDSlidingFitResult &slidingFitResult,
        pandora::CartesianVector &splitPosition) const;

    bool m_ignoreMichel;    ///< Whether to ignore transitions involving "michel" hits
    bool m_ignoreDiffuse;   ///< Whether to ignore transitions involving "diffuse" hits
};

} // namespace lar_content

#endif // LAR_SEMANTIC_SPLITTING_ALGORITHM_H
