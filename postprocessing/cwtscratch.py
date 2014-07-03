import matplotlib.pyplot as plt

spikeFilBinMean = np.mean(spikeFilBin, axis=0) / (binWidth*1e-3)
widths = np.round( np.arange(0.15,2.0,binWidth*dt/1000)/(binWidth*dt/1000) )
tmp = scipy.signal.find_peaks_cwt(spikeFilBinMean, widths, min_snr=0.5)
# low SNR was necessary to find "bad" peaks


plt.plot(bins[tmp]/1000.0, np.tile(1,np.shape(tmp)), 'gs')
plt.hold(True)
plt.plot( bins/1000.0, psthBin/(300*binWidth*1e-3), 'b-')
plt.plot( bins/1000.0, spikeFilBinMean, 'r-')
plt.show()
