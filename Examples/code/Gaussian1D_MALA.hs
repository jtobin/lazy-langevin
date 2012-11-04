import Numeric.MCMC.Langevin
import System.Random.MWC
import System.Environment
import System.Exit
import System.IO
import Numeric.AD
import Control.Monad

target :: RealFloat a => [a] -> a
target [m, v] = ((+ (- n * log v)) . (/ (-2 * v^2)) . sum . map ((^2) . (+ (- m)))) initData 
    where n        = fromIntegral (length initData)
          initData = [   3.4997144,   5.9977975,  -1.7228241,  -6.0734359, -10.7936751,  23.0753148
                     ,  19.0272481, -12.8223298,   6.3500001,   1.8686071, -15.7144080,   5.9247293
                     ,  -7.8529678,  11.3387758,   0.5488139,  -9.1889228,   4.6741005,  -4.5554079
                     ,  14.5353278,  -6.8592355,   7.2337063,   5.0775012,   3.8505277,  20.7062009
                     ,   3.1132428, -10.3894237, -11.5778168,   4.0464800,   7.9823569,  15.4726602 ]

gTarget :: [Double] -> [Double]
gTarget = grad target

main = do
    args  <- getArgs 
    when (args == []) $ do
        putStrLn  "(lazy-langevin) 1D Gaussian density                         "
        putStrLn  "Usage: ./Gaussian1D_MALA <numSteps> <stepSize> <inits>      " 
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

