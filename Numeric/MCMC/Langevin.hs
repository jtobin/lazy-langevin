{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.MCMC.Langevin where

import Control.Lens
import Control.Monad (when)
import Control.Monad.Primitive (PrimMonad, PrimState, RealWorld)
import Control.Monad.Trans.State.Strict (get, put, execStateT)
import qualified Data.Foldable as Foldable (foldl1)
import Data.Maybe (fromMaybe)
import Data.Sampling.Types
import Data.Traversable (for)
import Pipes hiding (for, next)
import qualified Pipes.Prelude as Pipes
import System.Random.MWC.Probability (Prob, Gen)
import qualified System.Random.MWC.Probability as MWC

-- | Trace 'n' iterations of a Markov chain and stream them to stdout.
--
-- >>> withSystemRandom . asGenIO $ mcmc 3 1 [0, 0] target
mcmc
  :: (Num (IxValue (t Double)), Show (t Double), Traversable t
     , FunctorWithIndex (Index (t Double)) t, Ixed (t Double)
     , IxValue (t Double) ~ Double)
  => Int
  -> Double
  -> t Double
  -> Target (t Double)
  -> Gen RealWorld
  -> IO ()
mcmc n step chainPosition chainTarget gen = runEffect $
        chain step Chain {..} gen
    >-> Pipes.take n
    >-> Pipes.mapM_ print
  where
    chainScore    = lTarget chainTarget chainPosition
    chainTunables = Nothing

-- A Markov chain driven by the Metropolis transition operator.
chain
  :: (Num (IxValue (t Double)), Traversable t
     , FunctorWithIndex (Index (t Double)) t, Ixed (t Double)
     , PrimMonad m, IxValue (t Double) ~ Double)
  => Double
  -> Chain (t Double) b
  -> Gen (PrimState m)
  -> Producer (Chain (t Double) b) m ()
chain step = loop where
  loop state prng = do
    next <- lift (MWC.sample (execStateT (langevin step) state) prng)
    yield next
    loop next prng

-- | Mean function for the discretized Langevin diffusion.
localMean
  :: (Fractional b, Ixed (f b), FunctorWithIndex (Index (f b)) f
     , IxValue (f b) ~ b)
  => Target (f b)
  -> b
  -> f b
  -> f b
localMean target e q = gzipWith (+) q scaled where
  scale  = (* (0.5 * e ^ (2 :: Int)))
  scaled = fmap scale (g q)
  g      = fromMaybe err (glTarget target)
  err    = error "adjustMomentum: no gradient provided"

perturb
  :: (Traversable f, PrimMonad m, Ixed (f Double)
     , FunctorWithIndex (Index (f Double)) f
     , IxValue (f Double) ~ Double)
  => Target (f Double)
  -> Double
  -> f Double
  -> Prob m (f Double)
perturb target e q = do
  zs  <- for q (const MWC.standard)
  let q0 = localMean target e q -- FIXME is there an error here?  looks dubious
      q1 = fmap (* e) zs
  return (gzipWith (+) q0 q1)

langevin
  :: (Floating (IxValue (t1 Double)), Ord (IxValue (t1 Double))
     , Traversable t1, PrimMonad m, Ixed (t1 Double)
     , FunctorWithIndex (Index (t1 Double)) t1
     , IxValue (t1 Double) ~ Double)
  => Double
  -> Transition m (Chain (t1 Double) b)
langevin e = do
  Chain {..} <- get
  proposal <- lift (perturb chainTarget e chainPosition)
  let currentMean   = localMean chainTarget e chainPosition
      proposalMean  = localMean chainTarget e proposal
      proposalScore = exp $ auxilliaryTarget chainTarget e
        (chainPosition, currentMean)
        (proposal, proposalMean)

      acceptProbability = whenNaN 0 proposalScore

  accept <- lift (MWC.bernoulli acceptProbability)
  when accept (put (Chain chainTarget proposalScore proposal chainTunables))

whenNaN :: RealFloat a => a -> a -> a
whenNaN val x
  | isNaN x   = val
  | otherwise = x

auxilliaryTarget
  :: (Floating (IxValue s), Ord (IxValue s), Foldable t, Foldable t1
     , Ixed s, FunctorWithIndex (Index s) t
     , FunctorWithIndex (Index s) t1, IxValue s ~ Double)
  => Target s
  -> IxValue s
  -> (s, t (IxValue s))
  -> (s, t1 (IxValue s))
  -> Double
auxilliaryTarget target e (q0, qm0) (q1, qm1) = min 0 $
    lTarget target q1 + log (isoGauss q1 qm0 (e ^ (2 :: Int)))
  - lTarget target q0 - log (isoGauss q0 qm1 (e ^ (2 :: Int)))

-- A container-generic zipwith.
gzipWith
  :: (FunctorWithIndex (Index s) f, Ixed s)
  => (a -> IxValue s -> b) -> f a -> s -> f b
gzipWith f xs ys = imap (\j x -> f x (fromMaybe err (ys ^? ix j))) xs where
  err = error "gzipWith: invalid index"

isoGauss
  :: (Floating (IxValue s), Ord (IxValue s), Foldable t, Ixed s
     , FunctorWithIndex (Index s) t)
  => s
  -> t (IxValue s)
  -> IxValue s
  -> IxValue s
isoGauss xs mu sig = exp (Foldable.foldl1 (+) (gzipWith density mu xs)) where
  density m x = log (densityNormal m sig x)
  densityNormal m s x
    | s < 0     = 0
    | otherwise =
          recip (sqrt (2 * pi) * s)
        * exp (negate (x - m) ^ (2 :: Int) / (2 * s ^ (2 :: Int)))

