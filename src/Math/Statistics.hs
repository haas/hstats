{-# OPTIONS_GHC -XBangPatterns #-}

-----------------------------------------------------------------------------
-- Module      : Math.Statistics
-- Copyright   : (c) 2008 Marshall Beddoe
-- License     : BSD3
--
-- Maintainer  : mbeddoe@<nospam>gmail.com
-- Stability   : experimental
-- Portability : portable
--
-- Description :
--   A collection of commonly used statistical functions.
-----------------------------------------------------------------------------

module Math.Statistics where

import Data.List
import Data.Ord (comparing)

-- Numerically stable mean
-- Thanks dmwit and ddarius for help on strictness issues
mean :: Floating a => [a] -> a
mean x = fst $ foldl' (\(!m, !n) x -> (m+(x-m)/(n+1),n+1)) (0,0) x

-- Harmonic mean
hmean :: (Floating a) => [a] -> a
hmean xs = fromIntegral (length xs) / (sum $ map (1/) xs)

-- Geometric mean
gmean :: (Floating a) => [a] -> a
gmean xs = (foldr1 (*) xs)**(1 / fromIntegral (length xs))

-- Median
median :: (Floating a, Ord a) => [a] -> a
median x | odd n  = head  $ drop (n `div` 2) x'
         | even n = mean $ take 2 $ drop i x'
                  where i = (length x' `div` 2) - 1
                        x' = sort x
                        n  = length x

-- Modes
-- Returns a sorted list of modes in descending order
modes :: (Ord a) => [a] -> [(Int, a)]
modes xs = sortBy (comparing $ negate.fst) $ map (\x->(length x, head x)) $ (group.sort) xs

-- Central moments
centralMoment xs 1 = 0
centralMoment xs r = (sum (map (\x -> (x-m)^r) xs)) / n
    where
      m = mean xs
      n = fromIntegral $ length xs

-- Range
range :: (Num a, Ord a) => [a] -> a
range xs = maximum xs - minimum xs

-- Average deviation
avgdev :: (Floating a) => [a] -> a
avgdev xs = mean $ map (\x -> abs(x - m)) xs
    where
      m = mean xs

-- Standard deviation
stddev :: (Floating a) => [a] -> a
stddev xs = sqrt $ var xs

-- Population variance
pvar :: (Floating a) => [a] -> a
pvar xs = centralMoment xs 2

-- Numerically stable sample variance
-- This crashes
var xs = (var' 0 0 0 xs) / (fromIntegral $ length xs - 1)
    where
      var' _ _ s [] = s
      var' m n s (x:xs) = var' nm (n + 1) (s + delta * (x - nm)) xs
         where
           delta = x - m
           nm = m + delta/(fromIntegral $ n + 1)

-- Interquartile range
-- XXX: Add case that takes into account even vs odd length
iqr xs = take (length xs - 2*q) $ drop q xs
    where
      q = ((length xs) + 1) `div` 4

-- Kurtosis
kurtosis xs = ((centralMoment xs 4) / (centralMoment xs 2)^2)-3

-- |Arbitrary quantile q of an unsorted list.  The quantile /q/ of /N/
-- |data points is the point whose (zero-based) index in the sorted
-- |data set is closest to /q(N-1)/.
quantile :: (Fractional b, Ord b) => Double -> [b] -> b
quantile q = quantileAsc q . sort

-- |As 'quantile' specialized for sorted data
quantileAsc :: (Fractional b, Ord b) => Double -> [b] -> b
quantileAsc _ [] = error "quantile on empty list"
quantileAsc q xs
    | q < 0 || q > 1 = error "quantile out of range"
    | otherwise = xs !! (quantIndex (length xs) q)
    where quantIndex :: Int -> Double -> Int
          quantIndex len q = case round $ q * (fromIntegral len - 1) of
                               idx | idx < 0    -> error "Quantile index too small"
                                   | idx >= len -> error "Quantile index too large"
                                   | otherwise  -> idx

-- Skew
skew xs = (centralMoment xs 3) / (centralMoment xs 2)**(3/2)

pearsonSkew1 xs = 3 * (mean xs - mo) / stddev xs
    where
      mo = snd $ head $ modes xs

pearsonSkew2 xs = 3 * (mean xs - median xs) / stddev xs

-- Covariance
cov :: (Floating a) => [a] -> [a] -> a
cov xs ys = sum (zipWith (*) (map f1 xs) (map f2 ys)) / (n - 1)
    where
      n = fromIntegral $ length $ xs
      m1 = mean xs
      m2 = mean ys
      f1 = \x -> (x - m1)
      f2 = \x -> (x - m2)

-- Covariance matrix
covMatrix :: (Floating a) => [[a]] -> [[a]]
covMatrix xs =  split' (length xs) cs
    where
      cs = [ cov a b | a <- xs, b <- xs]
      split' n = unfoldr (\y -> if null y then Nothing else Just $ splitAt n y)

-- Pearson's product-moment correlation coefficient
corr :: (Floating a) => [a] -> [a] -> a
corr x y = cov x y / (stddev x * stddev y)

-- |Least-squares linear regression of /y/ against /x/ for a
-- |collection of (/x/, /y/) data, in the form of (/b0/, /b1/, /r/)
-- |where the regression is /y/ = /b0/ + /b1/ * /x/ with Pearson
-- |coefficient /r/

linreg :: (Floating b) => [(b, b)] -> (b, b, b)
linreg xys = let !xs = map fst xys
                 !ys = map snd xys
                 !n = fromIntegral $ length xys
                 !sX = sum xs
                 !sY = sum ys
                 !sXX = sum $ map (^ 2) xs
                 !sXY = sum $ map (uncurry (*)) xys
                 !sYY = sum $ map (^ 2) ys
                 !alpha = (sY - beta * sX) / n
                 !beta = (n * sXY - sX * sY) / (n * sXX - sX * sX)
                 !r = (n * sXY - sX * sY) / (sqrt $ (n * sXX - sX^2) * (n * sYY - sY ^ 2))
             in (alpha, beta, r)
