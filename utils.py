import  matplotlib 
import  numpy              as      np
import  astropy.constants  as      const
import  pylab              as      pl


def set_size(w, h, ax=None):
    ''' 
    w, h: width, height in inches.
    '''
    if not ax: 
      ax = plt.gca()

    l    = ax.figure.subplotpars.left
    r    = ax.figure.subplotpars.right
    t    = ax.figure.subplotpars.top
    b    = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

def latexify(fig_width=None, fig_height=None, columns=1, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True):
    '''
    Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width:  float, optional, inches
    fig_height: float, optional, inches

    Columns : {1, 2}
    '''

    ##  Code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    ##  Width and max height in inches for IEEE journals taken from
    ##  computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    import  matplotlib.pyplot as      plt
    import  matplotlib        as      mpl

    from    numpy             import  sqrt

    assert(columns in [1,2])

    if ggplot:
      plt.style.use('ggplot')

    mpl.rc('text', usetex = usetex)

    if fig_width is None:
        fig_width = 6.08948 if columns == 2 else 6.08948 / 2.  ##  width in inches

    if fig_height is None:
      if equal is False:
        if ratio is None:
          ratio     = (sqrt(5.) - 1.0) / 2.0
        
        fig_height  = fig_width * ratio                        ## Height in inches

      else:
        fig_height  = fig_width 

    MAX_HEIGHT_INCHES = 12.1

    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + str(fig_height) + "so will reduce to" + str(MAX_HEIGHT_INCHES) + "inches.")

        fig_height = MAX_HEIGHT_INCHES

    params   = {'backend': 'pdf',
                ## 'text.latex.preamble': ['\usepackage{gensymb}'],
                'axes.labelsize': fontsize,  ## fontsize for x and y labels (was 10)
                'axes.titlesize': fontsize,
                'font.size': fontsize, 
                'legend.fontsize': fontsize,
                'xtick.labelsize': fontsize,
                'ytick.labelsize': fontsize,
                'text.usetex': usetex,
                'figure.figsize': [fig_width, fig_height],
                'font.family': 'serif',
                'figure.facecolor': 'w',
                'axes.facecolor': 'white'}

    matplotlib.rcParams.update(params)

    ax = pl.gca()

    ax.set_axis_on()

    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')

    if fig_width == None:
      w = 0.9 * 6.08948 if columns == 2 else 0.9 * 6.08948 / 2.

    else:
      w = fig_width 

    if equal is False:
      if ratio is None:
        ratio = (sqrt(5.) - 1.0) / 2.0

      h       = w * ratio                        ## Height in inches                                                                            

    else:
        h     =  w

    if fig_height != None:
      h       = fig_height 

    set_size(w, h, ax=ax)

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    '''
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.

    https://stackoverflow.com/questions/18311909/how-do-i-annotate-with-power-of-ten-formatting
    '''
    if not exponent:
        exponent = np.int(np.floor(np.log10(np.abs(num))))

    coeff = np.round(num / np.float(10.**exponent), decimal_digits)

    if not precision:
        precision = decimal_digits

    ##  \cdot    
    return r"${0:.{2}f} \times 10^{{{1:d}}}$".format(coeff, exponent, precision)


if __name__ == "__main__":
  print("\n\nWelcome to utils.\n\n")

  latexify(fig_width=None, fig_height=None, columns=1, equal=False)

  print("\n\nDone.\n\n")

