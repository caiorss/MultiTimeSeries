-- |
-- Module      :  MutliTimeSeries
-- Copyright   :  Robin S. Krom 2012
-- License     :  BSD3
-- 
-- Maintainer  :  Robin S. Krom
-- Stability   :  experimental
-- Portability :  unknown

module MultiTimeSeries (
        Sample
        , Vector
        , Matrix
        , WeightedSample
        , mean
        , stddev
        , crosscovariance
        , crossvariance
        , crosscorrelation
        , portmanteauTest
        , VarModel (..)
        , varModel
        , varForecast
        , varResidual
        , oslResidualCovariance
        , mlResidualCovariance
        , akaikeInfoCrit
        , weaklyStationary
        , mahalanobis
        , mahalanobisFilter
        , toSample
        , orderByTime
        , difference
)

where

import qualified Data.Packed as P
import qualified Numeric.Container as C
import qualified Numeric.LinearAlgebra as LA
import qualified Data.Vector.Generic as G
import qualified Data.Vector as V
import qualified Data.List.Zipper as Z
import Statistics.Distribution (Distribution, cumulative)
import Statistics.Distribution.ChiSquared
import Statistics.Test.Types (TestResult(..), significant)
import Data.List
import Data.Ord (comparing)
import Data.Complex (magnitude)

type Vector = P.Vector Double
type Matrix = P.Matrix Double
type Sample = [Vector]
type WeightedSample = [(Double, Vector)]
data VarModel = VarModel {phi0 :: Vector, phis :: [Matrix]} deriving (Show, Eq, Read) -- ^ (\phi_0, matrices \Phi_1 ... \Phi_p).


-- | Estimates the mean of a sample.
mean :: Sample -> Vector
mean s = C.scale (1 / fromIntegral (length s)) (sum s)

-- | Estimates the standart deviation.
stddev :: Sample -> Vector
stddev s = G.map sqrt $ C.scale (1 /fromIntegral (length s)) $ sum (map (G.map (^(2 :: Int)) . (\v -> v - mean s)) s)

-- | Estimates the crosscovariance matrix.
crosscovariance :: Sample -> Sample -> Matrix
crosscovariance s1 s2 = C.scale (1 / fromIntegral (length s1)) (sum $ zipWith C.outer s1 s2)

-- | Estimates the crossvariance with a certain lag of a sample
crossvariance :: Int -- ^ Lag.
                    -> Sample -- ^ The sample.
                    -> Matrix -- ^ Crossvariance matrix.
crossvariance l s = let 
    s' = map (\v -> v - mean s) s
    s'' = drop l s'
    in
    crosscovariance s' s''

-- | Estimates the crosscorrelation with a certain lag of a sample.
crosscorrelation :: Int -- ^ Lag. 
                    -> Sample -- ^ The sample.
                    -> Matrix -- ^ Crosscorrelation matrix.
-- Check for zero standart deviations.
crosscorrelation _ s | (C.prodElements . stddev) s == 0 = undefined
crosscorrelation l s = let
    m_D_inv = C.diag $ G.map (1/) $ stddev s
    m_Gamma_l = crossvariance l s
    in
    m_D_inv C.<> m_Gamma_l C.<> m_D_inv

-- | Multivariate Portmenteau test.
portmanteauTest :: Int -- ^ Dimension of the vector space. 
                    -> Int -- ^ Maximal lag.
                    -> Sample -- ^ The sample.
                    -> TestResult
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
            -> Sample -- ^ The sample. The size needs to be > p.
            -> VarModel -- ^ A var(p) model.
varModel dimV p s = let
    l = length s
    -- build matrices as described here:
    -- http://en.wikipedia.org/wiki/General_matrix_notation_of_a_VAR(p)
    m_Y = P.fromColumns $ drop p s
    m_Z = P.fromBlocks $ replicate (l - p) 1 : [map P.asColumn (take (l - p) (drop (p - i) s)) | i <- [1 .. p]]
    m_B = P.trans m_Z C.<\> P.trans m_Y

    cs_m_B = P.toColumns $ P.trans m_B

    submatrices [] = []
    submatrices cs = (take dimV cs : submatrices (drop dimV cs))

    in
    VarModel (head cs_m_B) (map P.fromColumns (submatrices $ tail cs_m_B))

-- | Forecast of a var(p)_model.
varForecast :: VarModel  -- ^ The var(p) model as returned by varModel.
                -> [Vector] -- ^ [r_{t-1}, r_{t-2}, ... r_{t-p}]
                -> Vector -- ^ The resulting forecast.
varForecast _ vs | null vs = undefined
varForecast vm rs = phi0 vm + sum (zipWith (C.<>) (phis vm) rs)

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
                        -> Sample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
oslResidualCovariance vm i s = C.scale (1/(fromIntegral (length s) - 2* fromIntegral i -1)) (residualCovariance vm i s)

-- | The residual covariance matrix estimate by maximum likelyhood method.
mlResidualCovariance :: VarModel
                        -> Int -- ^ The order of the residual
                        -> Sample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
mlResidualCovariance vm i s = C.scale (1/fromIntegral (length s)) (residualCovariance vm i s)

-- | Helper function for oslResidualCovariance and mlResidualCovariance.
residualCovariance :: VarModel -- ^ The var model.
                        -> Int -- ^ The order of the residual
                        -> Sample -- ^ Sample vec(r_1, r_2, r_3 ... ]
                        -> Matrix
residualCovariance vm i s = let
    l = length s
    a_i t = varResidual vm i (reverse (take t s))
    in
    sum [a_i t `C.outer` a_i t | t <- [i+1 .. l]]

-- | Select the var(p) model with the minimal akaike information criterion.
akaikeInfoCrit :: Int
                    -> Sample -- ^ The sample.
                    -> Int -- ^ Upper bound for the order of the var model.
                    -> VarModel -- ^ The best var model with minimal aic.
akaikeInfoCrit dimV s m = let
    aic vm = log (abs (LA.det $ mlResidualCovariance vm i s)) + (2.0 * fromIntegral dimV^(2 :: Int) * fromIntegral i) / fromIntegral (length s) where i = (length $ phis vm)
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
companionMatrix dimV vm = P.fromBlocks [[0, C.ident (dimV * (p - 1))], reverse $ phis vm] where p = length $ phis vm

type Norm = Vector -> Double

-- | Mahalanobis distance.
mahalanobis :: Matrix -- ^ Inverse covariance matrix.
               -> Norm
mahalanobis m_Sigma_inv v = v `C.dot` (m_Sigma_inv C.<> v)

distanceFilter :: Double -- ^ Scale factor
                    -> Double -- ^ Radius.
                    -> Norm  -- ^ Vector norm.
                    -> Sample -> Sample
distanceFilter s r n (v:vs) = v : go v vs
    where go x (y:xs) = if n (x - y) <= r 
                        then go x xs
                        else y : go (y + C.scale s (y - x)) xs
          go _ [] = []
distanceFilter _ _ _ [] = []

mahalanobisFilter :: Double  -- ^ Scale factor.
                        -> Double -- ^ Radius.
                        -> Matrix -- ^ Inverse covariance matrix.
                        -> Sample -> Sample
mahalanobisFilter s r m = distanceFilter s r (mahalanobis m)

-- | Calculate the trace of a matrix.
trace :: Matrix -> Double
trace m = G.sum $ LA.takeDiag m

-- | Create a sample from an array of 1-dimensional samples with timestamps.
toSample :: Ord a => V.Vector [(Double, a)] -> Sample
toSample = map (V.convert . V.map fst) . orderByTime

-- | Create a list of vectors with values and timestamps. Every list in the vector needs to contain at least one element.
orderByTime :: Ord b => V.Vector [(a, b)] -> [V.Vector (a, b)]
orderByTime xs = go zippers
    where 
            go zs
                    -- stop when we hit the end of a zipper.
                    | V.all (Z.endp . Z.right) zs = [V.map Z.cursor zs]
                    | otherwise = V.map Z.cursor zs : go (next zs)

            zippers = forwardTo maxi $ V.map Z.fromList xs
                where
                    maxi = V.maximum $ fmap (snd . head) xs

            forwardTo t zs
                        | V.all ((>= t) . snd . Z.cursor . Z.right) zs = zs
                        | otherwise = forwardTo t $ V.update zs $ V.map (\i -> (i, Z.right $ zs V.! i)) is
                            where is = V.findIndices ((< t) . snd . Z.cursor) zs

            next zs = V.update zs' $ V.map (\i -> (i, zs' V.! i)) minIndeces
                    where   
                        minIndeces = V.elemIndices mini zs''
                        -- zs' is not empty, otherwise go would have hit the first guard.
                        mini = V.minimum zs''
                        zs' = fmap safeRight zs
                        zs'' = fmap (snd . Z.cursor) zs'

            safeRight z = let z' = Z.right z in 
                if Z.endp z' then z else z'

-- | Take differences of a sample.
difference :: Sample -> Sample
difference (v : vs) = go v vs
    where   go x (y : ys) = (y - x) : go y ys
            go _ [] = []
difference [] = []
