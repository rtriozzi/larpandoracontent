/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/SemanticSplittingAlgorithm.h
 *
 *  @brief  Header file for the semantic splitting algorithm class.
 */

#ifndef LAR_SEMANTIC_SPLITTING_ALGORITHM_H
#define LAR_SEMANTIC_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

class SemanticSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    SemanticSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode DivideCaloHits(
        const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief Find the transition point for splitting the cluster based on hit semantic labels
     *
     *  @param pCluster the input cluster
     *  @param splitPosition the best layer
     */
    pandora::StatusCode FindBestSplitPosition(const pandora::Cluster *const pCluster, pandora::CartesianVector &splitPosition) const;
    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster,
                                       const float &splitPosition,
                                       pandora::CaloHitList &firstHitList,
                                       pandora::CaloHitList &secondHitList) const;

    // pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster,
    //                                    pandora::CaloHitList &firstHitList,
    //                                    pandora::CaloHitList &secondHitList) const override;

    bool m_ignoreMichel;    ///< Whether to ignore transitions involving Michel labels
    bool m_ignoreDiffuse;   ///< Whether to ignore transitions involving Diffuse labels
};

} // namespace lar_content

#endif // LAR_SEMANTIC_SPLITTING_ALGORITHM_H
