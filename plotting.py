import numpy as np
from matplotlib import pyplot as plt
from pyshellspec.auxilliary import add_wing_phase


def plot_light_curve(outputname, obsphase, obsmag, error, synmag, whole_lc_phase, whole_lc_mag, passband=''):
    """
    Plots a light curve comparison.
    :param outputname:
    :param obsphase:
    :param obsmag:
    :param error:
    :param synmag:
    :param whole_lc_phase:
    :param whole_lc_mag:
    :param passband:
    :return:
    """
    # pad the synthetic light curve
    pad_lc_phase, pad_lc_mag = add_wing_phase(whole_lc_phase, whole_lc_mag)

    # First plot comparison of the LC observed and synthetic
    fig = plt.figure(figsize=(15, 10))
    ax = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    ax.errorbar(obsphase, obsmag, yerr=error, color='r', fmt='^')
    ax.plot(pad_lc_phase, pad_lc_mag, 'k-')
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(1.05 * max([obsmag.max(), whole_lc_mag.max()]), 0.95 * min([obsmag.min(), whole_lc_mag.min()]))
    ax.set_xlabel('Orbital phase')
    ax.set_ylabel('%s(mag)' % passband)

    # compute residuals
    residuals = obsmag - synmag

    # Second plot the residuals
    ax = plt.subplot2grid((4, 1), (3, 0), rowspan=1)
    ax.errorbar(obsphase, residuals, yerr=error, color='y', fmt='o')
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(1.05 * residuals.min(), 1.05 * residuals.max())
    ax.set_xlabel('Orbital phase')
    ax.set_ylabel('residuals(mag)')

    # plot and save it
    fig.tight_layout()
    fig.savefig(outputname)
    plt.close()


def plot_squared_visibility(figname, u, v, vis2data, vis2err, vis2syn, image=None, scale=None):
    """
    Plots comparison of observed and synthetic squared visibility.
    :param figname:
    :param u:
    :param v:
    :param vis2data:
    :param vis2err:
    :param vis2syn:
    :return:
    """
    # compute baseline
    baseline = np.sqrt(u ** 2 + v ** 2)

    # divide data based on position angle
    posang = np.degrees(np.arctan(v / u))
    idns = np.where(np.abs(posang) > 60.)
    idew = np.where(np.abs(posang) < 30.)
    ides = np.where((np.abs(posang) > 30.) & (np.abs(posang) < 60.))

    # do the plotting -- first the observations and synthetics
    fig = plt.figure(figsize=(15, 12))
    ax = plt.subplot2grid((6, 5), (0, 0), rowspan=3, colspan=5)
    ax.errorbar(baseline[idew], vis2data[idew], color='r', fmt='^', yerr=vis2err[idew], label='PA < 30 deg')
    ax.errorbar(baseline[ides], vis2data[ides], color='g', fmt='^', yerr=vis2err[ides], label='PA > 30 deg,PA < 60 deg')
    ax.errorbar(baseline[idns], vis2data[idns], color='b', fmt='^', yerr=vis2err[idns], label='PA > 60 deg')
    ax.plot(baseline, vis2syn, 'ko', label='Synthetic')
    ax.set_xlabel('Spatial frequency(rad$^{-1}$)')
    ax.set_ylabel('$V^2$')
    ax.set_ylim(-0.1, 1.2)
    # plt.legend()

    # second the residuals
    residual = vis2data - vis2syn

    ax = plt.subplot2grid((6, 5), (3, 0), rowspan=1, colspan=5)
    ax.errorbar(baseline[idew], residual[idew], yerr=vis2err[idew], color='r', fmt='o')
    ax.errorbar(baseline[ides], residual[ides], yerr=vis2err[ides], color='g', fmt='o')
    ax.errorbar(baseline[idns], residual[idns], yerr=vis2err[idns], color='b', fmt='o')
    ax.set_xlabel('Spatial frequency(rad$^{-1}$)')
    ax.set_ylabel('$residuals$')
    ax.set_ylim(1.05 * residual.min(), 1.05 * residual.max())

    # plot the uv coverage
    ax = plt.subplot2grid((6, 5), (4, 0), rowspan=2, colspan=2)
    ax.plot(u[idew], v[idew], 'r.')
    ax.plot(u[ides], v[ides], 'g.')
    ax.plot(u[idns], v[idns], 'b.')
    ax.plot(-u[idew], -v[idew], 'r.')
    ax.plot(-u[ides], -v[ides], 'g.')
    ax.plot(-u[idns], -v[idns], 'b.')

    # set the limit
    lim = max([np.abs(u).max(), np.abs(v).max()])
    ax.set_xlim(-1.05 * lim, 1.05 * lim)
    ax.set_ylim(-1.05 * lim, 1.05 * lim)

    ax.set_xlabel('u(rad$^{-1}$)')
    ax.set_ylabel('v(rad$^{-1}$)')
    ax.set_aspect('equal')
    # print u, v

    if image is not None:
        ax = plt.subplot2grid((6, 5), (4, 3), rowspan=2, colspan=2)
        ax.imshow(image, cmap='gray', extent=[scale.min(), scale.max(), scale.min(), scale.max()])
        ax.set_xlim(scale.min(), scale.max())
        ax.set_ylim(scale.min(), scale.max())
        ax.set_xlabel(r'$\alpha$(rad)')
        ax.set_ylabel(r'$\delta$(rad)')

    # save the image
    fig.tight_layout()
    fig.savefig(figname)
    plt.close()


def plot_triple_product(figname, uc1, vc1, uc2, vc2, t3amp, t3phi, t3amperr, t3phierr, t3ampsyn, t3phisyn):
    """
    Plots comparison between observed and synthetic triple product.
    :param figname:
    :param uc1:
    :param vc1:
    :param uc2:
    :param vc2:
    :param t3amp:
    :param t3phi:
    :param t3amperr:
    :param t3phierr:
    :param t3ampsyn:
    :param t3phisyn:
    :return:
    """
    # compute the third baselines
    uc3 = uc1 + uc2
    vc3 = vc1 + vc2

    # triple baseline
    t3bas = np.sqrt(uc1 ** 2 + vc1 ** 2 + uc2 ** 2 + vc2 ** 2 + uc3 ** 2 + vc3 ** 2)

    # plot the triple amplitude comparison
    fig = plt.figure(figsize=(15, 10))
    ax = plt.subplot2grid((4, 4), (0, 0), rowspan=3, colspan=2)
    ax.errorbar(t3bas, t3amp, yerr=t3amperr, color='r', fmt='^', label=figname)
    ax.plot(t3bas, t3ampsyn, 'ko', label=figname)
    ax.set_xlabel('Spatial frequency (rad$^{-1}$)')
    ax.set_ylabel(r'$\vert T_3\vert$')
    ax.set_ylim(-0.2, 1.6)

    # plot the triple phase comparison
    ax = plt.subplot2grid((4, 4), (0, 2), rowspan=3, colspan=2)
    ax.errorbar(t3bas, t3phi, yerr=t3phierr, color='r', fmt='^', label=figname)
    ax.plot(t3bas, t3phisyn, 'ko', label=figname)
    ax.set_xlabel('Spatial frequency (rad$^{-1}$)')
    ax.set_ylabel(r'$T_3\phi$(deg)')
    ax.set_ylim(1.05 * min([t3phi.min(), t3phisyn.min()]),
                1.05 * max([t3phi.max(), t3phisyn.max()]))

    # plot triple amplitude residuals
    residual = t3amp - t3ampsyn
    ax = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=2)
    ax.errorbar(t3bas, residual,  yerr=t3amperr, color='y', fmt='o', label=figname)
    ax.set_xlabel('Spatial frequency (rad$^{-1}$)')
    ax.set_ylabel(r'$residuals$')
    ax.set_ylim(-1.6, 1.6)

    # plot triple phase residuals
    residual = t3phi - t3phisyn
    ax = plt.subplot2grid((4, 4), (3, 2), rowspan=1, colspan=2)
    ax.errorbar(t3bas, residual,  yerr=t3phierr, color='y', fmt='o', label=figname)
    ax.set_xlabel('Spatial frequency (rad$^{-1}$)')
    ax.set_ylabel(r'$residuals$(deg)')
    ax.set_ylim(1.05 * min([residual.min()]), 1.05 * max([residual.max()]))

    # save the image
    fig.tight_layout()
    fig.savefig(figname)
    plt.close()
