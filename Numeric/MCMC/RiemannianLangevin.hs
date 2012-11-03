{-# OPTIONS_GHC -Wall #-}

module Numeric.MCMC.RiemannianLangevin (
            MarkovChain(..)
          , Parameters, createParameters
          , runChain
          ) where

import Numeric.MCMC.Langevin hiding (Parameters, runChain)
import Numeric.LinearAlgebra  
import Control.Arrow
import Control.Monad
import Control.Monad.Primitive
import Control.Monad.Reader
import System.Random.MWC
import System.Random.MWC.Distributions 

-- | Parameters for the chain.  The target (expected to be a log density), the
--   right matrix product of its curvature and gradient, the inverse Fisher metric tensor,
--   its (Cholesky) square root, and the step size.
data Parameters = Parameters { 
    _target              :: [Double] -> Double        -- Target (log density)
  , _curvatureXgradient  :: [Double] -> [Double]      -- Curvature right-multiplied by gradient
  , _invFisherMetric     :: [Double] -> Matrix Double -- Inverse Fisher metric tensor
  , _sqrtInvFisherMetric :: [Double] -> Matrix Double -- Cholesky component of tensor
  , _eps                 :: {-# UNPACK #-} !Double    -- Step size
  }

type ViewsParameters = ReaderT Parameters

-- | Construct Parameters (data constructor not exported).
createParameters :: ([Double] -> Double)     -- Target (log density)
                 -> ([Double] -> [Double])   -- Gradient
                 -> ([Double] -> [[Double]]) -- Hessian
                 -> Double                   -- Step size
                 -> Parameters               -- Parameters
createParameters t g h = 
    Parameters t curvatureXgradient invFisherMetric sqrtInvFisherMetric 
  where curvatureXgradient xs = 
            let mat = invFisherMetric xs <> fromColumns [fromList (g xs)]
            in  concat . toLists . trans $ mat
        invFisherMetric     = inv . fromRows . map (fromList . map (* (-1))) . h
        sqrtInvFisherMetric = trans . chol . invFisherMetric
{-# INLINE createParameters #-}

-- | Non-isotropic Gaussian density.
nonIsoGauss :: [Double] -> [Double] -> Matrix Double -> Double
nonIsoGauss xs mu sig = exp val
  where val = -0.5*p*log (2*pi) - 0.5*ldet - 0.5*
            (trans (xsM `sub` muM) <> invSig <> (xsM `sub` muM)) @@> (0, 0)
        (xsM, muM) = (\f (a, b) -> (f a, f b)) (\l -> fromColumns [fromList l]) (xs, mu)
        p                   = fromIntegral $ cols sig
        (invSig, (ldet, _)) = invlndet sig

-- | Mean function for the discretized Riemannian Langevin diffusion.
localMean :: [Double]               -- Current state
          -> ([Double] -> [Double]) -- Curvature * gradient
          -> Double                 -- Step size
          -> [Double]               -- Localized mean of proposal distribution
localMean t c e = zipWith (+) t (map (* (0.5 * e^(2 :: Int))) (c t))
{-# INLINE localMean #-}

-- | Perturb the state, creating a new proposal.
perturb :: PrimMonad m 
        => [Double]                      -- Current state
        -> Gen (PrimState m)             -- MWC PRNG
        -> ViewsParameters m [Double]    -- Resulting perturbation.
perturb t g = do
    Parameters _ c _ s e <- ask
    zs <- replicateM (length t) (lift $ standard g)
    let adjustedBrownianMotion = s t <> fromColumns [fromList zs]
        abmList = concat . toLists . trans $ adjustedBrownianMotion
        perturbedState = zipWith (+) t0 t1 
        t0 = localMean t c e
        t1 = map (* eps) abmList
    return $! perturbedState 
{-# INLINE perturb #-}

-- | Perform a Metropolis accept/reject step.
metropolisStep :: PrimMonad m 
               => MarkovChain                   -- Current state
               -> Gen (PrimState m)             -- MWC PRNG
               -> ViewsParameters m MarkovChain -- New state
metropolisStep state g = do
    Parameters target c iF _ e <- ask
    let (t0, nacc) = (theta &&& accepts) state
    zc       <- lift $ uniformR (0, 1) g
    proposal <- perturb t0 g
    let mc = if   zc < acceptProb 
             then (proposal, 1)
             else (t0,       0)

        acceptProb = if isNaN val then 1 else val where val = arRatio 
        
        arRatio = exp . min 0 $ target proposal - target t0 
          + log (nonIsoGauss t0 (localMean proposal c e) ((e^(2 :: Int)) `scale` iF proposal))
          - log (nonIsoGauss proposal (localMean t0 c e) ((e^(2 :: Int)) `scale` iF t0))

    return $! Config (fst mc) (nacc + snd mc)
{-# INLINE metropolisStep #-}

-- | Diffuse through states.
runChain :: Parameters      -- Parameters of the Markov chain.
         -> Int             -- Number of epochs to iterate the chain.
         -> MarkovChain     -- Initial state of the Markov chain.
         -> Gen RealWorld   -- MWC PRNG
         -> IO MarkovChain  -- End state of the Markov chain, wrapped in IO.
runChain params nepochs initConfig g 
    | nepochs == 0 = return initConfig
    | otherwise    = do
        result <- runReaderT (metropolisStep initConfig g) params
        print result
        runChain params (nepochs - 1) result g

