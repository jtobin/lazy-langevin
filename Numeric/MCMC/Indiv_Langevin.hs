import Numeric.AD
import System.IO (stderr, hPutStrLn)
import System.Environment (getArgs)
import Control.Error
import Control.Monad
import Control.Monad.Primitive
import Control.Arrow hiding (left)
import System.Random.MWC
import System.Random.MWC.Distributions
import Data.List
import Statistics.Distribution
import Statistics.Distribution.Normal hiding (standard)
import Control.Proxy
import Control.Monad.Trans
import qualified Test.QuickCheck as QC

-- "globals"
target [x0, x1] = (-1)*(5*(x1 - x0^2)^2 + 0.2*(1 - x0)^2)
target _        = error "explode"
{-# INLINE target #-}

gradTarget :: [Double] -> [Double]
gradTarget = grad target
{-# INLINE gradTarget #-}

eps = 0.755

tt = [1.0, 1.5] :: [Double]
tp = [0.9778541492302911,1.4499536814731855] :: [Double]

-- funcs

-- | State of the Markov Chain.  Current parameter values are held in 'theta', while
--   accepts counts the number of proposals accepted.
data MarkovChain = Config { theta   :: [Double] 
                          , accepts :: {-# UNPACK #-} !Int }

-- | Display the current state (but not 
instance Show MarkovChain where
    show config = filter (`notElem` "[]") $ show (theta config)

-- | Density function for an isotropic Gaussian.  The (identity) covariance 
--   matrix is multiplied by the scalar 'sig'.
isotropicGaussian :: [Double] -> [Double] -> Double -> Double
isotropicGaussian xs mu sig = foldl1' (*) (zipWith density nds xs)
    where nds = map (`normalDistr` sig) mu
{-# INLINE isotropicGaussian #-}

-- | Perturb the state, creating a new proposal.
-- NOTE tested
perturb :: PrimMonad m => [Double] -> Gen (PrimState m) -> m [Double]
perturb theta g = do
    zs <- replicateM (length theta) (standard g)
    return $! perturbedState theta zs
  where perturbedState :: [Double] -> [Double] -> [Double]
        perturbedState theta zs = zipWith (+) t0 t1
            where t0 = localizedProposalMean theta
                  t1 = map (* eps) zs
{-# INLINE perturb #-}

-- | The mean function for the discretized Langevin diffusion.
-- NOTE move into scope of acceptProb
-- NOTE depends on eps, target
localizedProposalMean ::  [Double] -> [Double]
localizedProposalMean theta = zipWith (+) theta (map (* (0.5 * eps^2)) (gradTarget theta))
{-# INLINE localizedProposalMean #-}

-- | Calculate the accept/reject ratio
-- NOTE move into scope of acceptProb
-- NOTE depends on eps, target
arRatio :: [Double] -> [Double] -> Double
arRatio theta proposal = exp . min 0 $
    target proposal + log (isotropicGaussian theta (localizedProposalMean proposal) (eps^2))
  - target theta    - log (isotropicGaussian proposal (localizedProposalMean theta) (eps^2))
{-# INLINE arRatio #-}

-- | The acceptance probability of a proposal.
-- NOTE possibly move into scope of metropolisComparison
acceptProb :: [Double] -> [Double] -> Double
acceptProb theta proposal | isNaN val = 1
                          | otherwise = val
    where val = arRatio theta proposal 
{-# INLINE acceptProb #-}

prop_AcceptProbNonNeg :: [Double] -> [Double] -> Bool
prop_AcceptProbNonNeg xs ys = acceptProb xs ys >= 0

-- | Accept or reject a proposal, given a comparison probability.
-- NOTE possibly move into scope of metropolisStep
metropolisComparison :: [Double] -> [Double] -> Double -> ([Double], Int)
metropolisComparison theta proposal zc = if   zc < acceptProb theta proposal
                                         then (proposal, 1)
                                         else (theta,    0)
{-# INLINE metropolisComparison #-}

-- | A metropolis accept/reject step.
metropolisStep :: PrimMonad m => MarkovChain -> Gen (PrimState m) -> m MarkovChain
metropolisStep state g = do
    let (t0, nacc) = (theta &&& accepts) state
    zc       <- uniformR (0, 1) g
    proposal <- perturb t0 g
    let mc = metropolisComparison t0 proposal zc
    return $! Config (fst mc) (nacc + snd mc)

-- | Diffuse through states.
runChain :: Int -> MarkovChain -> Gen RealWorld -> IO MarkovChain
runChain nepochs initConfig g | nepochs == 0 = do
    hPutStrLn stderr ("(" ++ show (accepts initConfig) ++ " accepts)")
    return initConfig
                              | otherwise    = do
    result <- metropolisStep initConfig g
    print result
    runChain (nepochs - 1) result g

main :: IO ()
main = do
    [n0, e0] <- getArgs 
    let (nepochs, eps) = (read n0, read e0) :: (Int, Double)
    g <- create
    runChain nepochs (Config [0, 0] 0) g
    return ()






-- chain :: Gen RealWorld -> MarkovChain -> Producer MarkovChain IO r
-- chain g t0 = forever $ do
--     z <- lift $ metropolisStep t0 g
--     yield z
-- 
-- taker :: Int -> Pipe a a IO ()
-- taker n = replicateM_ n $ do
--     z <- await
--     yield z
-- 
-- printer :: Show a => Consumer a IO r
-- printer = forever $ do
--     z <- await
--     lift $ print z
-- 
-- testPipeline = do 
--     g <- create
--     let t0 = Config [0.0, 0.0] 0 
--     runPipe $ printer <+< taker 10 <+< chain g t0
--
-- main = do
--     g <- create
--     runChain (Config [0.0, 0.0] 0) g
--
-- rnds :: Gen RealWorld -> Producer Double IO r
-- rnds g = forever $ do
--     z <- lift $ uniformR (0, 1) g
--     yield z
--
