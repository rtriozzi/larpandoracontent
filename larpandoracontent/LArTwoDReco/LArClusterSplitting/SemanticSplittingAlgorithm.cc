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

StatusCode SemanticSplittingAlgorithm::DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstHitList, pandora::CaloHitList &secondHitList) const
{
    float splitPositionX(0.f);
    firstHitList.clear();
    secondHitList.clear();

    if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(pCluster, splitPositionX)) {
        StatusCode status = this->DivideCaloHits(pCluster, splitPositionX, firstHitList, secondHitList);
        std::cout << "Semantic split applied at x = " << splitPositionX 
                  << " | hits: " << firstHitList.size() << " + " << secondHitList.size() << std::endl;
        return status;
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SemanticSplittingAlgorithm::FindBestSplitPosition(
    const pandora::Cluster *const pCluster, float &splitPositionX) const
{

    std::vector<std::tuple<CartesianVector, std::string, float>> hitInfo; ///< store position, label, confidence
    const OrderedCaloHitList &orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for (auto iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector hits(iter->second->begin(), iter->second->end());
        for (const CaloHit *const pCaloHit : hits)
        {
            const auto *pLArHit = static_cast<const lar_content::LArCaloHit *>(pCaloHit);
            const auto &allScores = pLArHit->GetHitScores();
            const auto &allLabels = pLArHit->GetHitScoreLabels();

            if (allScores.empty() || allLabels.empty())
                continue;

            std::vector scores(allScores.begin() + 1, allScores.end());
            std::vector labels(allLabels.begin() + 1, allLabels.end());

            auto sortedScores = scores;
            std::sort(sortedScores.begin(), sortedScores.end(), std::greater<float>());
            const float confidence = sortedScores[0] / sortedScores[1];
            // if (confidence < 1.5f)
            //     continue;

            const size_t bestIdx = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));

            hitInfo.emplace_back(pLArHit->GetPositionVector(), labels[bestIdx], confidence);
        }
    }

    if (hitInfo.size() < 2)
        return STATUS_CODE_NOT_FOUND;    

    std::sort(hitInfo.begin(), hitInfo.end(),
              [](const auto &a, const auto &b) { return std::get<0>(a).GetX() < std::get<0>(b).GetX(); });

    std::vector<float> transitionPoints;
    std::vector<float> transitionScores;

    for (size_t i = 1; i < hitInfo.size(); ++i)
    {
        const std::string &prevLabel = std::get<1>(hitInfo[i - 1]);
        const std::string &currLabel = std::get<1>(hitInfo[i]);
        if ((m_ignoreMichel && (prevLabel == "michel" || currLabel == "michel")) ||
            (m_ignoreDiffuse && (prevLabel == "diffuse" || currLabel == "diffuse")))
            continue;

        if (prevLabel != currLabel)
        {
            const float splitX = 0.5f * (std::get<0>(hitInfo[i - 1]).GetX() + std::get<0>(hitInfo[i]).GetX());
            transitionPoints.push_back(splitX);

            const float avgConf = 0.5f * (std::get<2>(hitInfo[i - 1]) + std::get<2>(hitInfo[i]));
            transitionScores.push_back(avgConf);

            std::cout << "Transition between labels " << prevLabel << " and " << currLabel
                        << " at x = " << splitX << " with average confidence " << avgConf << std::endl;
        }
    }

    std::cout << "Found " << transitionPoints.size() << " candidate splitting points." << std::endl;

    if (transitionPoints.empty())
        return STATUS_CODE_NOT_FOUND;

    const size_t bestIdx = std::distance(
        transitionScores.begin(), std::max_element(transitionScores.begin(), transitionScores.end()));

    splitPositionX = transitionPoints[bestIdx];

    return STATUS_CODE_SUCCESS;

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SemanticSplittingAlgorithm::DivideCaloHits(const pandora::Cluster *const pCluster, const float &splitPositionX, pandora::CaloHitList &firstHitList, pandora::CaloHitList &secondHitList) const
{
    const OrderedCaloHitList &orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for (auto iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector hits(iter->second->begin(), iter->second->end());
        for (const CaloHit *const pCaloHit : hits)
        {
            if (pCaloHit->GetPositionVector().GetX() < splitPositionX)
                firstHitList.push_back(pCaloHit);
            else
                secondHitList.push_back(pCaloHit);
        }
    }

    std::cout << "Cluster split into two lists with #hits: " << firstHitList.size() << " and " << secondHitList.size() << std::endl;
    
    if (firstHitList.empty() || secondHitList.empty())
        return STATUS_CODE_NOT_FOUND;
    
    if (firstHitList.size() < 3 || secondHitList.size() < 3)
        return STATUS_CODE_NOT_FOUND;

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
