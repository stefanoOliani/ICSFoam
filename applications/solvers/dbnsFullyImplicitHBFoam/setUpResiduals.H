residualsIO initResI(2,1);

initResI.sInitRes()[0] = initRes().sInitRes()[2*K];
initResI.sInitRes()[1] = initRes().sInitRes()[2*K+1];
initResI.vInitRes()[0] = initRes().vInitRes()[K];

initResI.sFinalRes()[0] = initRes().sFinalRes()[2*K];
initResI.sFinalRes()[1] = initRes().sFinalRes()[2*K+1];
initResI.vFinalRes()[0] = initRes().vFinalRes()[K];

solnControl[K].setResidual(initResI);
