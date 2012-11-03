import Numeric.MCMC.Langevin
import System.Random.MWC
import System.Environment
import System.IO
import Numeric.AD

target :: RealFloat a => [a] -> a
target [x0, x1] = (-1)*(5*(x1 - x0^2)^2 + 0.05*(1 - x0)^2)

gradTarget :: [Double] -> [Double]
gradTarget = grad target

main = do
    [n0, e0] <- getArgs
    let nepochs = read n0 :: Int
        eps     = read e0 :: Double
    g <- create
    let params  = Parameters target gradTarget eps
        config0 = Config [0.0, 0.0] 0
    results <- runChain params nepochs config0 g
    hPutStrLn stderr $ "(" ++ show (accepts results) ++ " accepts)"
    return ()

