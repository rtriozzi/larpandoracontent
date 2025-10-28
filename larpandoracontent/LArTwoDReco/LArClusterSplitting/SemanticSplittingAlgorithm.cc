/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/SemanticSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the semantic splitting algorithm class.
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/SemanticSplittingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

SemanticSplittingAlgorithm::SemanticSplittingAlgorithm() :
    m_ignoreMichel(true),
    m_ignoreDiffuse(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode SemanticSplittingAlgorithm::FindBestSplitPosition(
    const TwoDSlidingFitResult &slidingFitResult, pandora::CartesianVector &splitPosition) const
{
    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    struct HitInfo {
        pandora::CartesianVector position;  ///< hit position
        std::string label;                  ///< assigned semantic label
        float confidence;                   ///< confidence in the assigned label
        float rL;                           ///< longitudinal position along the sliding fit direction
    };
    std::vector<HitInfo> hits;

    for (auto iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (const CaloHit *pCaloHit : *(iter->second))
        {
            const auto *pLArCaloHit = static_cast<const lar_content::LArCaloHit *>(pCaloHit);
            const auto &allScores = pLArCaloHit->GetHitScores();
            const auto &allLabels = pLArCaloHit->GetHitScoreLabels();

            if (allScores.empty() || allLabels.empty())
                continue;

            std::vector<float> scores(allScores.begin() + 1, allScores.end());
            std::vector<std::string> labels(allLabels.begin() + 1, allLabels.end());

            // Skip hits with low confidence in their best predicted label
            std::vector<float> sortedScores = scores;
            std::sort(sortedScores.begin(), sortedScores.end(), std::greater<float>());
            const float confidence = sortedScores[0] / sortedScores[1];
            if (confidence < 1.5f)
                continue;

            // Get the best label for this hit
            const size_t bestIdx = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));

            // Project onto sliding fit direction
            float rL(0.f), rT(0.f);
            slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), rL, rT);

            hits.push_back({pCaloHit->GetPositionVector(), labels[bestIdx], confidence, rL});
        }
    }

    if (hits.size() < 2)
        return STATUS_CODE_NOT_FOUND;

    std::sort(hits.begin(), hits.end(), 
        [](const HitInfo &a, const HitInfo &b) { return a.rL < b.rL; });
    std::vector<float> transitionConfidence;
    std::vector<size_t> transitionIndices;

    for (size_t i = 1; i < hits.size(); ++i)
    {
        const std::string &prevLabel = hits[i - 1].label;
        const std::string &currLabel = hits[i].label;

        // Skip diffuse and Michel hits, if requested
        if ((m_ignoreMichel && (prevLabel == "michel" || currLabel == "michel")) ||
            (m_ignoreDiffuse && (prevLabel == "diffuse" || currLabel == "diffuse")))
            continue;

        // Store transition points between different labels
        if (prevLabel != currLabel)
        {
            transitionConfidence.push_back(0.5f * (hits[i - 1].confidence + hits[i].confidence));
            transitionIndices.push_back(i);
        }
    }

    if (transitionConfidence.empty())
        return STATUS_CODE_NOT_FOUND;

    // Split at the transition with the highest average confidence
    const size_t bestIdx = std::distance(
        transitionConfidence.begin(), std::max_element(transitionConfidence.begin(), transitionConfidence.end()));
    splitPosition = 
        (hits[transitionIndices[bestIdx] - 1].position 
        + hits[transitionIndices[bestIdx]].position)
        * 0.5f;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SemanticSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IgnoreMichel", m_ignoreMichel));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IgnoreDiffuse", m_ignoreDiffuse));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
