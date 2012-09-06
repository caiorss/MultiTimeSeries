-- |
-- Module      :  MutliTimeSeries
-- Copyright   :  Robin S. Krom 2012
-- License     :  BSD3
-- 
-- Maintainer  :  Robin S. Krom
-- Stability   :  experimental
-- Portability :  unknown


module MultiTimeSeries (
MSample,
Vector,
Matrix,
MWeightedSample,
mean,
crossvariance,
crosscorrelation,
portmanteauTest,
varModel,
varForecast,
varResidual,
oslResidualCovariance,
mlResidualCovariance,
akaikeInfoCrit,
weaklyStationary,
mahalanobisDistance,
filterMahalanobis,
prop_mean_invariant_under_reversion,
prop_mean_linear,
prop_cv_symmetric,
prop_cv_transformation,
prop_ccv_switch_args,
prop_ccv_sum
)

where

import qualified Data.Packed as P
import qualified Numeric.Container as C
import qualified Numeric.LinearAlgebra as LA
import qualified Data.Vector.Generic as G
import Statistics.Distribution (Distribution, cumulative)
import Statistics.Distribution.ChiSquared
import Statistics.Test.Types (TestResult(..), significant)
import Data.List (minimumBy)
import Data.Ord (comparing)
import Data.Complex (magnitude)

type Vector = P.Vector Double
type Matrix = P.Matrix Double
type MSample = [Vector]
type MWeightedSample = [(Double, Vector)]
type VarModel = (Vector, [Matrix]) -- ^ (\phi_0, matrices \Phi_1 ... \Phi_p).

-- | Rounding error upper bound. Can get terrible if matrices are of bad shape in tests.
roundingError :: Double 
roundingError = 1e-7

-- | Estimates the mean of a sample.
mean :: MSample -> Vector
mean s | badMSample s = undefined
mean s = C.scale (1 / fromIntegral (length s)) (sum s)

prop_mean_invariant_under_reversion :: MSample -> Bool
prop_mean_invariant_under_reversion s = C.norm1 (mean s - mean (reverse s)) < roundingError

prop_mean_linear :: (Double, MSample, MSample) -> Bool
prop_mean_linear (mu, s1, s2) = C.norm1 (C.scale mu (mean s1) + mean s2 - mean (zipWith (+) (map (C.scale mu) s1) s2)) < roundingError

-- | Estimates the crosscovariance matrix.
crosscovariance :: MSample -> MSample -> Matrix
crosscovariance s1 s2 | badMSample s1 || badMSample s2 = undefined
crosscovariance s1 s2 = C.scale (1 / fromIntegral (length s1)) (sum $ zipWith C.outer s1 s2)

prop_ccv_switch_args :: (MSample, MSample) -> Bool
prop_ccv_switch_args (s1, s2) = C.norm1 (P.flatten (crosscovariance s1 s2 - P.trans (crosscovariance s2 s1))) < roundingError

-- terrible rounding errors.
prop_ccv_sum :: (MSample, MSample) -> Bool
prop_ccv_sum (s1, s2) = C.norm1 (P.flatten (crosscovariance s s - crosscovariance s1 s1 - crosscovariance s1 s2 - crosscovariance s2 s1 - crosscovariance s2 s2)) < 10 where s = zipWith (+) s1 s2

-- | Estimates the crossvariance with a certain lag of a sample
crossvariance :: Int -- ^ Lag.
                    -> MSample -- ^ The sample.
                    -> Matrix -- ^ Crossvariance matrix.
crossvariance _ s | badMSample s = undefined
crossvariance l s = let 
    s' = map (\v -> v - mean s) s
    s'' = drop l s'
    in
    crosscovariance s' s''

prop_cv_symmetric :: MSample -> Bool
prop_cv_symmetric s = C.norm1 (P.flatten (m - P.trans m)) < roundingError where m = crossvariance 0 s

-- terrible rounding errors.
prop_cv_transformation :: (MSample, Matrix, Vector) -> Bool
prop_cv_transformation (s, m_A, a) = C.norm1 (P.flatten (crossvariance 0 (map (\x -> m_A C.<> x + a) s) - (m_A C.<> crossvariance 0 s C.<> P.trans m_A))) < 10

-- | Estimates the crosscorrelation with a certain lag of a sample.
crosscorrelation :: Int -- ^ Lag. 
                    -> MSample -- ^ The sample.
                    -> Matrix -- ^ Crosscorrelation matrix.
crosscorrelation _ s | badMSample s = undefined
-- Check for zero standart deviation.
crosscorrelation _ s | C.prodElements (sum (map (G.map (^(2 :: Int)) . (\v -> v - mean s)) s)) == 0 = undefined
crosscorrelation l s = let
    s' = map (\v -> v - mean s) s
    m_D_inv = C.diag $ G.map ((1/) . sqrt) (sum $ map (G.map (^(2 :: Int))) s')
    m_Gamma_l = crossvariance l s
    in
    m_D_inv C.<> m_Gamma_l C.<> m_D_inv

-- | Multivariate Portmenteau test.
portmanteauTest :: Int -- ^ Dimension of the vector space. 
                    -> Int -- ^ Maximal lag.
                    -> MSample -- ^ The sample.
                    -> TestResult
portmanteauTest _ _ s | badMSample s = undefined
portmanteauTest dimV _ _ | dimV == 0 = undefined
portmanteauTest dimV m s = let
    m_Gamma_0_inv = LA.inv $ crosscorrelation 0 s 
    m_Gamma i = crosscorrelation i s
    s_T = length s
    q m' = sum [(fromIntegral s_T^(2 :: Int) / fromIntegral (s_T - l)) * trace (C.trans (m_Gamma l) C.<> m_Gamma_0_inv C.<> m_Gamma l C.<> m_Gamma_0_inv) | l <- [1 .. m']]
    in
    significant $ cumulative (chiSquared (dimV * m^(2 :: Int))) (q m) < 0.05

-- | Estimate the var(p) model of a given sample according to
-- r_t = \phi_0 + \sum_{i=1}^p \Phi_i r_{t-i} + e_t
varModel :: Int -- ^ Dimension of the vector space.
            -> Int -- ^ p parameter of the model.
            -> MSample -- ^ The sample.
            -> VarModel -- ^ A var(p) model.
varModel _ _ s | badMSample s = undefined
varModel dimV p s = let
    l = length s
    -- build matrices as described here:
    -- http://en.wikipedia.org/wiki/General_matrix_notation_of_a_VAR(p)
    m_Y = P.fromColumns $ drop (l - p) s
    m_Z = P.fromBlocks ([1] : [map P.asColumn (take (l - p) (drop (p - i) s)) | i <- [1 .. p]])
    m_B = m_Y C.<\> m_Z
    cs_m_B = P.toColumns m_B
    submatrices cs = (take dimV cs : submatrices (drop dimV cs))
    in
    (head cs_m_B, map P.fromColumns (submatrices $ tail cs_m_B))

-- | Forecast of a var(p)_model.
varForecast :: VarModel  -- ^ The var(p) model as returned by varModel.
                -> [Vector] -- ^ [r_{t-1}, r_{t-2}, ... r_{t-p}]
                -> Vector -- ^ The resulting forecast.
varForecast _ vs | null vs = undefined
varForecast vm rs = fst vm + sum (zipWith (C.<>) (snd vm) rs)

-- | The i-th residual of a var(p) model.
varResidual :: VarModel -- ^ The var(p) model as returned by varModel.
                -> Int -- ^ Order of the residual.
                -> [Vector] -- ^ [r_t, r_{t-1}, r_{t-2}, ... r_{t-i}]
                -> Vector -- ^ The i-th residual, a^i_t.
varResidual _ _ vs | null vs = undefined
varResidual vm  i rs = head rs - varForecast vm (take i (tail rs))

-- | The residual covariance matrix estimate by osl method.
oslResidualCovariance :: VarModel
                        -> Int -- ^ The order of the residual
                        -> MSample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
oslResidualCovariance _ _ s | badMSample s = undefined
oslResidualCovariance vm i s = C.scale (1/(fromIntegral (length s) - 2* fromIntegral i -1)) (residualCovariance vm i s)

-- | The residual covariance matrix estimate by maximum likelyhood method.
mlResidualCovariance :: VarModel
                        -> Int -- ^ The order of the residual
                        -> MSample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
mlResidualCovariance _ _ s | badMSample s = undefined
mlResidualCovariance vm i s = C.scale (1/fromIntegral (length s)) (residualCovariance vm i s)

-- | Helper function for oslResidualCovariance and mlResidualCovariance.
residualCovariance :: VarModel -- ^ The var model.
                        -> Int -- ^ The order of the residual
                        -> MSample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
residualCovariance _ _ s | badMSample s = undefined
residualCovariance vm i s = let
    l = length s
    a_i t = varResidual vm i (reverse (take t s))
    in
    sum [a_i t `C.outer` a_i t | t <- [i+1 .. l]]

-- | Select the var(p) model with the minimal akaike information criterion.
akaikeInfoCrit :: Int
                    -> MSample -- ^ The sample.
                    -> Int -- ^ Upper bound for the order of the var model.
                    -> VarModel -- ^ The best var model with minimal aic.
akaikeInfoCrit _ s _ | badMSample s = undefined
akaikeInfoCrit dimV s m = let
    aic vm = log (abs (LA.det $ mlResidualCovariance vm i s)) + (2.0 * fromIntegral dimV^(2 :: Int) * fromIntegral i) / fromIntegral (length s) where i = (length $ snd vm)
    in
    minimumBy (comparing aic) [varModel dimV i s | i <- [0 .. m]]

-- | Determines wether a given var(p) model is weakly stationary by looking 
-- the biggest eigenvalue of the companion matrix.
weaklyStationary :: Int -- ^ Dimension of vector space.
                    -> VarModel -- ^ The var model.
                    -> Bool
weaklyStationary i _ | i <= 0 = undefined
weaklyStationary dimV vm = maxEV < 1 where maxEV = G.maximum $ G.map magnitude (LA.eigenvalues $ companionMatrix dimV vm)

-- | Builds the companion matrix for a var(p) model.
companionMatrix :: Int -- ^ Dimension of vector space.
                    -> VarModel -- ^ The var model.
                    -> Matrix
companionMatrix i _ | i <= 0 = undefined
companionMatrix dimV vm = P.fromBlocks [[0, C.ident (dimV * (p - 1))], reverse $ snd vm] where p = length $ snd vm

-- | Mahalanobis distance.
mahalanobisDistance :: Matrix -- ^ Inverse covariance matrix.
                        -> Vector
                        -> Double
mahalanobisDistance m_Sigma_inv v = v `C.dot` (m_Sigma_inv C.<> v)

-- | Filter by Mahalanobis distance.
filterMahalanobis :: Matrix -- ^ Inverse covariance matrix.
                        -> Double -- ^ Radius.
                        -> [Vector] -> [Maybe Vector] -- ^ Filter function.
filterMahalanobis m_Sigma mu vs = h vs Nothing where
    h [] _ = []
    h (y:ys) Nothing = Just y : h ys (Just y)
    h (y:ys) (Just c) = if mahalanobisDistance m_Sigma (y - c) < mu then Nothing : h ys (Just c) else Just y : h ys (Just y)

-- | Calculate the trace of a matrix.
trace :: Matrix -> Double
trace m = G.sum $ LA.takeDiag m 

badMSample :: MSample -> Bool
badMSample s = null s || G.null (head s)
