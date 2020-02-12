import matplotlib
import numpy             as np
import matplotlib.pyplot as plt

from   mpl_toolkits.axes_grid1 import   make_axes_locatable


def fast_scatter(ax, xs, ys, values, mmin, mmax, N, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0):
  ##  Digitize.
  step         = (mmax - mmin) / N
  
  points       =  mmin + step * np.floor((np.clip(values, a_min=mmin, a_max=mmax) - mmin) / step)
  points      -=  1.e-1

  levels       = np.unique(points)
  
  if printit:
    print(mmin, mmax, N, step)
    print(np.unique(values))
    print(levels)

    print(len(xs), len(ys), len(values))
    
  ##
  indexs       = np.floor((np.clip(levels, a_min=mmin, a_max=mmax) - mmin) / step).astype(np.int)
    
  if len(levels) > 500:
    raise  ValueError('{}'.format(levels))

  ##  ['Summer', 'Wistia', 'autumn_r']
  cmap         = plt.get_cmap(cmap, N)
  norm         = matplotlib.colors.Normalize(vmin=mmin, vmax=mmax)

  colors       = cmap([1. * x / N for x in np.arange(0, N, 1)])

  for i, level in enumerate(levels):
    isin       = (points == level)
    
    if printit:
      print('Plotting level {} of {} - {} targets at {}.'.format(i, len(levels), np.count_nonzero(isin), level))

    ax.plot(xs[isin], ys[isin], markersize=markersize, c=colors[indexs[i]], lw=0, marker='.', alpha=alpha)

  ##
  ##  divider    = make_axes_locatable(ax)

  ##  cax        = divider.append_axes('right', size='2%', pad=0.2)                                                                                                                                                 
  ##  cb         = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)

  return  0
