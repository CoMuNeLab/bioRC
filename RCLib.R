library(compiler)

################################################
# Setup the Reservoir and run it
################################################
runReservoir.f <- function(Signal, trainLen, testLen, initLen, W, inSize=1, outSize=1, RC_leakRate = 0.3, RC_mode = "GENERATIVE", RC_reg = 1e-8, RC_bias_intensity=0.1){
    #Signal             Input signal
    #trainLen           Number of samples for training
    #testLen            Number of samples for testing
    #initLen            
    #W                  Weighted adjacency matrix of the reservoir
    #inSize             Number of input signals
    #outSize            Number of output signals
    
    #RC_leakRate        leaking rate
    #RC_mode            GENERATIVE | PREDICTIVE
    #RC_reg             regularization coefficient

    RC_Nodes <- dim(W)[1]

    #############################
    # generate the ESN reservoir
    #############################

    Win = matrix(runif(RC_Nodes*(1+inSize),-0.5,0.5),RC_Nodes)

    # normalizing and setting spectral radius (correct, slow):
    cat('Computing spectral radius...')
    rhoW = abs(eigen(W,only.values=TRUE)$values[1])
    cat('done.\n')
    W = W * 1.25 / rhoW

    # allocated memory for the design (collected states) matrix
    X = matrix(0,1+inSize+RC_Nodes,trainLen-initLen)
    # set the corresponding target matrix directly
    Yt = matrix(data[(initLen+2):(trainLen+1)],1)

    bias <- runif(RC_Nodes,-0.5,0.5)*RC_bias_intensity #runif(n, min = 0, max = 1)

    #############################
    # run the reservoir with the data and collect X
    #############################
    x = rep(0,RC_Nodes)
    for (t in 1:trainLen){
      u = data[t]
      x = (1-RC_leakRate)*x + RC_leakRate*tanh( Win %*% rbind(1,u) + W %*% x + bias )
      if (t > initLen)
        X[,t-initLen] = rbind(1,u,x)
    }

    #############################
    # train the output
    #############################
    X_T = t(X)
    Wout = Yt %*% X_T %*% solve( X %*% X_T + RC_reg*diag(1+inSize+RC_Nodes) )

    # run the trained ESN in a generative mode. no need to initialize here, 
    # because x is initialized with training data and we continue from there.
    Y = matrix(0,outSize,testLen)
    u = data[trainLen+1]

    if(RC_mode=="PREDICTIVE"){
      for (t in 1:testLen){
        x = (1-RC_leakRate)*x + RC_leakRate*tanh( Win %*% rbind(1,u) + W %*% x + bias )
        y = Wout %*% rbind(1,u,x)
        Y[,t] = y

        u = data[trainLen+t+1] 
      }
    }else{
      for (t in 1:testLen){
        x = (1-RC_leakRate)*x + RC_leakRate*tanh( Win %*% rbind(1,u) + W %*% x + bias )
        y = Wout %*% rbind(1,u,x)
        Y[,t] = y

        u = y
      }  
    }

    return(list(Signal=data, Output=Y, W=W, Activations=X, W.in=Win, W.out=Wout))
    #return(list(Signal=data, Output=Y, W.out=Wout))
}
runReservoir <- cmpfun(runReservoir.f)

################################################
# Compute MSE for the first errorLen time steps
################################################
getMSE <- function(RC, trainLen, errorLen, norm=F){
    #RC                 Obtained from runReservoir
    #trainLen           Number of samples for training
    #errorLen           Number of samples in testing for estimating the error

    mse = ( sum( (RC$Signal[(trainLen+2):(trainLen+errorLen+1)] - RC$Output[1,1:errorLen])^2 )
            / errorLen )

    if(norm){
        mse = mse/sd(RC$Signal[(trainLen+2):(trainLen+errorLen+1)])
    }

    return(mse)
}

addWhiteNoise <- function(data, SNRdb){
    return( as.matrix(data + runif(nrow(data),-1,1)*sqrt(3)*sd(data)*10^(-SNRdb/20)) )  
}

SNRdb <- function(data, noise){
    return( 10*log(sd(data)^2/sd(noise)^2,10) )
}

normSignal <- function(data){
    return( (data-min(data))/(max(data)-min(data)) - 0.5 )
}


