%% test case 1: uniform spectrum and 

testSphere  = ita_generateSampling_equiangular(10,10);
 
testAudioObj = itaAudio(ones(256,testSphere.nPoints),44100,'freq');
testAudioObj.channelCoordinates = testSphere;

[DFcompSphere,commSphere] = ita_diffuseFieldCompensation(testAudioObj);

testHrtfObj = itaHRTF(merge(testAudioObj,testAudioObj));
[DFcompSphereHrtf,commSphereHrtf] = ita_diffuseFieldCompensation(testHrtfObj);

assert(all(abs(commSphere.freqData-1)<eps,'all'),'Common component should be equal to unity')
assert(all(abs(commSphereHrtf.freqData-1)<eps,'all'),'Common component should be equal to unity')

assert(all(abs(DFcompSphere.freqData-1)<eps,'all'),'DTF component should be equal to unity')
assert(all(abs(DFcompSphereHrtf.freqData-1)<eps,'all'),'DTF component should be equal to unity')

%%  test case 2 -> itaHRTF from ITA artificial head
dummyHead = ita_read('\\verdi\fileserver\database\measurements\HRTFs\[HRTF measurements]\openSource\finishedHRTF_5deg.ita');
[DFcompDummy,commDummy] = ita_diffuseFieldCompensation(dummyHead);


%% validate against Ramonas old implementation (now only found in user scripts)

[DFcompDummy_rbo,commDummy_rbo] = test_rbo_DTF_itaHRTF(dummyHead);


%% check differences (differences <.05dB to allow for compitational inaccuracies)
assert(max(abs(commDummy_rbo.ch(1:2).freqData_dB-commDummy.ch(1:2).freqData_dB),[],'all')<0.05,'Results from old RBO implementation and new should be the same for equiangular sampling')

%% plot difference
difference = ita_divide_spk(commDummy,commDummy_rbo.ch(1:2));
difference.pf
ylim([-1,1])
