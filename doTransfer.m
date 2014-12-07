tstart = tic;
T1_tf = deformationTransfer3(S0, T0, S1);
ttrans = toc(tstart);
fprintf('transferred in %f seconds.\n', ttrans);
showComp;