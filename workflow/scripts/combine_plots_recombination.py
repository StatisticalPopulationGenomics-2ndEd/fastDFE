"""
Combine plots of SFS stratified by rho and rho histogram.
"""
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1.1, 1]})  # 2 rows, 1 column

img1 = mpimg.imread('sfs.rho.png')
img2 = mpimg.imread('rho.hist.png')

axs[0].imshow(img1)
axs[0].axis('off')

axs[1].imshow(img2)
axs[1].axis('off')

fig.text(0.06, 0.96, 'A', fontsize=16, fontweight='bold')
fig.text(0.06, 0.48, 'B', fontsize=16, fontweight='bold')

plt.tight_layout(pad=0)
plt.savefig('sfs.hist.rho.png', dpi=500)

plt.show()
