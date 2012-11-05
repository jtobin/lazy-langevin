{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE BangPatterns #-}

module Numeric.MCMC.Langevin ( 
          MarkovChain(..), Options(..)
        , runChain
        ) where

import Control.Monad
import Control.Monad.Trans
import Control.Monad.Reader
import Control.Monad.Primitive
import Control.Arrow 
import System.Random.MWC
import System.Random.MWC.Distributions
import Data.List
import Statistics.Distribution
import Statistics.Distribution.Normal hiding (standard)
import GHC.Float

-- | State of the Markov chain.  Current parameter values are held in 'theta', 
--   while accepts counts the number of proposals accepted.
data MarkovChain = MarkovChain { theta   :: [Double] 
                               , accepts :: {-# UNPACK #-} !Int }

-- | Options for the chain.  The target (expected to be a log density), its 
--   gradient, and a step size tuning parameter.
data Options = Options { _target  :: [Double] -> Double
                       , _gTarget :: [Double] -> [Double]
                       , _eps     :: {-# UNPACK #-} !Double }

-- | A result with this type has a view of the chain options.
type ViewsOptions = ReaderT Options

-- | Display the current state. 
instance Show MarkovChain where
    show config = filter (`notElem` "[]") $ show (map double2Float (theta config))

-- | Density function for an isotropic Gaussian.  The (identity) covariance 
--   matrix is multiplied by the scalar 'sig'.
isoGauss :: [Double] -> [Double] -> Double -> Double
isoGauss xs mu sig = foldl1' (*) (zipWith density nds xs)
    where nds = map (`normalDistr` sig) mu
{-# INLINE isoGauss #-}

-- | Mean function for the discretized Langevin diffusion.
localMean :: Monad m 
          => [Double]                -- Current state
          -> ViewsOptions m [Double] -- Localized mean 
localMean t = do
    Options _ gTarget e <- ask
    return $! zipWith (+) t (map (* (0.5 * e^(2 :: Int))) (gTarget t))
{-# INLINE localMean #-}

-- | Perturb the state, creating a new proposal.
perturb :: PrimMonad m 
        => [Double]                   -- Current state
        -> Gen (PrimState m)          -- MWC PRNG
        -> ViewsOptions m [Double]    -- Resulting perturbation.
perturb t g = do
    Options _ _ e <- ask
    zs    <- replicateM (length t) (lift $ standard g)
    t0    <- localMean t
    let perturbedState = zipWith (+) t0 t1 
        t1 = map (* e) zs
    return $! perturbedState 
{-# INLINE perturb #-}

-- | Perform a Metropolis accept/reject step.
metropolisStep :: PrimMonad m 
               => MarkovChain                -- Current state
               -> Gen (PrimState m)          -- MWC PRNG 
               -> ViewsOptions m MarkovChain -- New state
metropolisStep state g = do
    Options target _ e <- ask
    let (t0, nacc) = (theta &&& accepts) state
    zc           <- lift $ uniformR (0, 1) g
    proposal     <- perturb t0 g
    t0Mean       <- localMean t0
    proposalMean <- localMean proposal
    let mc = if   zc < acceptProb 
             then (proposal, 1)
             else (t0,       0)

        acceptProb = if isNaN val then 0 else val where val = arRatio 
        
        arRatio = exp . min 0 $ 
            target proposal + log (isoGauss proposal t0Mean (e^(2 :: Int)))
          - target t0       - log (isoGauss t0 proposalMean (e^(2 :: Int)))

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

