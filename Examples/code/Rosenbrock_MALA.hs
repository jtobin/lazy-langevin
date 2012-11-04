import Numeric.MCMC.Langevin
import System.Random.MWC
import System.Environment
import System.Exit
import System.IO
import Numeric.AD
import Control.Monad

target :: RealFloat a => [a] -> a
target [x0, x1] = (-1)*(5*(x1 - x0^2)^2 + 0.05*(1 - x0)^2)

gTarget :: [Double] -> [Double]
gTarget = grad target

main = do
    args  <- getArgs 
    when (args == []) $ do
        putStrLn  "(lazy-langevin) Rosenbrock density                          "
        putStrLn  "Usage: ./Rosenbrock_MALA <numSteps> <stepSize> <inits>      " 
        putStrLn  "                                                            "
        putStrLn  "numSteps         : Number of Markov chain iterations to run."
        putStrLn  "stepSize         : Perturbation scaling parameter.          "
        putStrLn  "inits            : Filepath containing points at which to   "
        putStrLn  "                   initialize the chain.                    "
        exitSuccess

    inits <- fmap (map read . words) (readFile (args !! 2)) :: IO [Double]

    let nepochs = read (head args) :: Int
        eps     = read (args !! 1) :: Double
        params  = Options target gTarget eps
        config  = MarkovChain inits 0

    g       <- create
    results <- runChain params nepochs 1 config g

    hPutStrLn stderr $ 
        let nAcc  = accepts results
            total = nepochs 
        in  show nAcc ++ " / " ++ show total ++ " (" ++ 
              show ((fromIntegral nAcc / fromIntegral total) :: Float) ++ 
              ") proposals accepted"

