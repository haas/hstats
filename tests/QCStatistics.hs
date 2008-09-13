module QCStatistics where

import Control.Monad
import Test.QuickCheck

import Math.Statistics

-- Quantiles

prop_Quantiles :: Property
prop_Quantiles = forAll (choose (0, 1)) $
                 \q -> forAll (sized $ \n -> replicateM (n + 1) (arbitrary :: Gen Double)) $
                 \xs -> let nx = fromIntegral $ length xs
                            qx = quantile q xs
                        in collect (count (< qx) xs, count (<= qx) xs, count (> qx) xs, count (>= qx) xs) $
                           and [count (< qx) xs <= ceiling (q * nx),
                                count (<= qx) xs >= floor (q * nx),
                                count (> qx) xs <= ceiling ((1 - q) * nx),
                                count (>= qx) xs >= floor ((1 - q) * nx)]
    where count pred = length . filter pred

-- Linear regression

prop_LinReg :: (Double, Double) -> Property
prop_LinReg (a0, a1) = forAll genXYs $
                       \xys -> let (b0, b1, r) = linreg xys
                               in and [b0 ~= a0, b1 ~= a1,
                                       if a1 ~= 0.0 then isNaN r else abs r ~= 1]
    where genXYs :: Gen [(Double, Double)]
          genXYs = liftM (map (\x -> (x, a0 + a1 * x))) $ genXs
          genXs :: Gen [Double]
          genXs = liftM (scanl1 (+)) $ sized $ \n -> replicateM (n + 2) $ genNonZero
          genNonZero :: Gen Double
          genNonZero = ap (elements [id,negate]) $ choose (1, 10) 

(~=) :: (Fractional a, Ord a) => a -> a -> Bool
x1 ~= x2 = abs(x1 - x2) < 1e-6             
