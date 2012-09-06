{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

import Test.QuickCheck
import qualified Data.Vector.Generic as G
import qualified Data.Packed as P
import MultiTimeSeries 


arbitraryMSample :: Int -> Gen MSample
arbitraryMSample i = listOf1 (fmap G.fromList (vectorOf i arbitrary))

arbitraryMSample' :: Int -> Int -> Gen MSample
arbitraryMSample' l d = vectorOf l (fmap G.fromList (vectorOf d arbitrary))

arbitraryVector :: Int -> Gen Vector
arbitraryVector i = fmap G.fromList (vectorOf i arbitrary)

arbitrarySquareMatrix :: Int -> Gen Matrix
arbitrarySquareMatrix i = fmap P.fromLists (vectorOf i (vectorOf i arbitrary))

newtype MSPair = MSP (MSample, MSample) deriving (Show)

unwrapMSP :: MSPair -> (MSample, MSample)
unwrapMSP (MSP (s1, s2)) = (s1, s2)

instance Arbitrary MSPair where
    arbitrary = do
        l <- elements [1 .. 10]
        d <- elements [1 .. 5]
        s1 <- arbitraryMSample' l d
        s2 <- arbitraryMSample' l d
        return $ MSP (s1, s2)

samples :: Gen MSample
samples = do
        l <- elements [1 .. 5]
        arbitraryMSample l
    
main :: IO ()
main = do
    putStrLn "mean invariant under reversion."
    quickCheck $ forAll samples prop_mean_invariant_under_reversion
    putStrLn "\nmean linear."
    quickCheck (prop_mean_linear . (\(mu, MSP (s1, s2)) -> (mu, s1, s2)))
    putStrLn "\ncrossvariance with lag 0 symmetric."
    quickCheck $ forAll samples prop_cv_symmetric
    putStrLn "\ncrosscovariance transform correct when switching arguments."
    quickCheck (prop_ccv_switch_args . unwrapMSP)
    putStrLn "\ncrosscovariance of a sum of samples."
    quickCheck (prop_ccv_sum . unwrapMSP)
    putStrLn "\ncrossvariance transforms correct."
    quickCheck (forAll triple prop_cv_transformation)
        where triple = do
                        i <- elements [1 .. 5]
                        s <- arbitraryMSample i
                        m <- arbitrarySquareMatrix i
                        v <- arbitraryVector i
                        return (s,m,v)
