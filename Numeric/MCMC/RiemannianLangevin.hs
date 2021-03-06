{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE BangPatterns #-}

module Numeric.MCMC.RiemannianLangevin (
            MarkovChain(..)
          , Options, createOptions
          , runChain, localMean, perturb
          ) where

import Numeric.MCMC.Langevin hiding (Options, runChain)
import Numeric.LinearAlgebra  
import Control.Arrow
import Control.Monad
import Control.Monad.Primitive
import Control.Monad.Reader
import System.Random.MWC
import System.Random.MWC.Distributions 

-- | Options for the chain.  The target (expected to be a log density), the
--   right matrix product of its curvature and gradient, the inverse Fisher metric tensor,
--   its (Cholesky) square root, and the step size.
data Options = Options { 
    _target              :: [Double] -> Double        -- Target (log density)
  , _curvatureXgradient  :: [Double] -> [Double]      -- Curvature right-multiplied by gradient
  , _invFisherMetric     :: [Double] -> Matrix Double -- Inverse Fisher metric tensor
  , _sqrtInvFisherMetric :: [Double] -> Matrix Double -- Square root of the tensor
  , _eps                 :: {-# UNPACK #-} !Double    -- Step size
  }

-- | A result with this type has a view of the chain options.
type ViewsOptions = ReaderT Options

-- | Construct Options (data constructor not exported).
-- FIXME sqrtm is not converging; alternatives?
createOptions :: ([Double] -> Double)     -- Target (log density)
              -> ([Double] -> [Double])   -- Gradient
              -> ([Double] -> [[Double]]) -- Hessian
              -> Double                   -- Step size
              -> Options                  -- Options
createOptions t g h = 
    Options t curvatureXgradient invFisherMetric sqrtInvFisherMetric 
  where curvatureXgradient xs = 
            let mat = invFisherMetric xs <> fromColumns [fromList (g xs)]
            in  concat . toLists . trans $ mat
        invFisherMetric       = pinv . fromRows . map (fromList . map (* (-1))) . h 
        sqrtInvFisherMetric x = let z = schur (invFisherMetric x) 
                                in  fst z <> sqrtm (snd z) <> trans (fst z)
{-# INLINE createOptions #-}

-- | Non-isotropic Gaussian density.
nonIsoGauss :: [Double] -> [Double] -> Matrix Double -> Double
nonIsoGauss xs mu sig = exp val
  where val = -0.5*p*log (2*pi) - 0.5*ldet - 0.5*
            (trans (xsM `sub` muM) <> invSig <> (xsM `sub` muM)) @@> (0, 0)
        (xsM, muM) = (\f (a, b) -> (f a, f b)) (\l -> fromColumns [fromList l]) (xs, mu)
        p                   = fromIntegral $ cols sig
        (invSig, (ldet, _)) = invlndet sig
{-# INLINE nonIsoGauss #-}

-- | Mean function for the discretized Riemannian Langevin diffusion.
localMean :: Monad m 
          => [Double]                -- Current state
          -> ViewsOptions m [Double] -- Localized mean of proposal distribution
localMean t = do
    Options _ c _ _ e <- ask
    return $! zipWith (+) t (map (* (0.5 * e^(2 :: Int))) (c t))
{-# INLINE localMean #-}

-- | Perturb the state, creating a new proposal.
perturb :: PrimMonad m 
        => [Double]                   -- Current state
        -> Gen (PrimState m)          -- MWC PRNG
        -> ViewsOptions m [Double]    -- Resulting perturbation.
perturb t g = do
    Options _ _ _ s _ <- ask
    zs <- replicateM (length t) (lift $ standard g)
    t0 <- localMean t
    let adjustedBrownianMotion = s t <> fromColumns [fromList zs]
        abmList = concat . toLists . trans $ adjustedBrownianMotion
        perturbedState = zipWith (+) t0 t1 
        t1 = map (* eps) abmList
    return $! perturbedState 
{-# INLINE perturb #-}

-- | Perform a Metropolis accept/reject step.
metropolisStep :: PrimMonad m 
               => MarkovChain                -- Current state
               -> Gen (PrimState m)          -- MWC PRNG
               -> ViewsOptions m MarkovChain -- New state
metropolisStep state g = do
    Options target _ iF _ e <- ask
    let (t0, nacc) = (theta &&& accepts) state
    zc           <- lift $ uniformR (0, 1) g
    proposal     <- perturb t0 g
    t0Mean       <- localMean t0
    proposalMean <- localMean proposal
    let mc = if   zc < acceptProb 
             then (proposal, 1)
             else (t0,       0)

        acceptProb = if isNaN val then 1 else val where val = arRatio 
        
        arRatio = exp . min 0 $ 
            target proposal + log (nonIsoGauss proposal t0Mean ((e^(2::Int)) `scale` iF t0)) 
          - target t0       - log (nonIsoGauss t0 proposalMean ((e^(2::Int)) `scale` iF proposal)) 

    return $! MarkovChain (fst mc) (nacc + snd mc)
{-# INLINE metropolisStep #-}

-- | Diffuse through states.
runChain :: Options         -- Options of the Markov chain.
         -> Int             -- Number of epochs to iterate the chain.
         -> Int             -- Print every nth iteration
         -> MarkovChain     -- Initial state of the Markov chain.
         -> Gen RealWorld   -- MWC PRNG
         -> IO MarkovChain  -- End state of the Markov chain, wrapped in IO.
runChain = go
  where go o n t !c g | n == 0 = return c
                      | n `rem` t /= 0 = do
                            r <- runReaderT (metropolisStep c g) o
                            go o (n - 1) t r g
                      | otherwise = do
                            r <- runReaderT (metropolisStep c g) o
                            print r
                            go o (n - 1) t r g
{-# INLINE runChain #-}

