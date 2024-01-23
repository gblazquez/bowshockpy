import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib

from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import matplotlib.font_manager as fm

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from itertools import product

from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.coordinates import SkyCoord

from photutils.isophote import EllipseGeometry
from photutils import EllipticalAperture, EllipticalAnnulus, CircularAperture, RectangularAperture




def mosaic_slices(image_cube, nrow, ncol, chan_0, chan_f, wcs, box,
                  operation=None, width_chan=0, output_name=None,
                  ellip_dic=None, return_ims=False, points_bluearm=None,
                  **kwargs):
    """
    Makes a mosaic of spectral slices of the cube
    box = [[x_1,y_1],[x_2,y_2]]
    """

    for param in mf.default_params:
        kwargs[param] = kwargs[param] \
                if param in kwargs \
                else mf.default_params[param]
#   changes the size of labels
    matplotlib.rcParams.update({'font.size': kwargs['font_size']})
#   set False to unable latex rendering
    plt.rc('text', usetex=kwargs['use_tex'])
#   plt.rc('font', family='serif')

    if kwargs['offset_coordinates']:
        center_skycoords = SkyCoord(mf.default_params["vla4b_deg"][0]*u.deg,
                                    mf.default_params["vla4b_deg"][1]*u.deg)

        center_pixels = [
            float(skycoord_to_pixel(center_skycoords, wcs)[0]),
            float(skycoord_to_pixel(center_skycoords, wcs)[1])
        ]

        x0y0_extent = mf.pixel2offset([0,0],
                             [kwargs["header"]["CDELT1"]*3600,
                              kwargs["header"]["CDELT2"]*3600],
                             center_pixels,
                             bounds=True)

        xNyN_extent = mf.pixel2offset([kwargs["header"]["NAXIS1"]-1,
                                        kwargs["header"]["NAXIS2"]-1],
                             [kwargs["header"]["CDELT1"]*3600,
                              kwargs["header"]["CDELT2"]*3600],
                             center_pixels,
                             bounds=True)

        extent = [x0y0_extent[0], xNyN_extent[0],
                   x0y0_extent[1], xNyN_extent[1]]
        x0y0lim = mf.pixel2offset(box[0],
                             [kwargs['header']["CDELT1"]*3600,
                              kwargs['header']["CDELT2"]*3600],
                             center_pixels,
                             bounds=True)

        xNyNlim = mf.pixel2offset(box[1],
                             [kwargs['header']["CDELT1"]*3600,
                              kwargs['header']["CDELT2"]*3600],
                             center_pixels,
                             bounds=True)

    else:
        extent = None

    if kwargs['subgridspec'] is not None:
        supgrid = kwargs['subgridspec']
        gs = supgrid.subgridspec(
            nrow, ncol+2,
            wspace=kwargs['wspace'],
            hspace=kwargs['hspace'],
            height_ratios=[1]*nrow,
            width_ratios=[1]*ncol + [kwargs["blankratio"], kwargs["axcbar_ratio"]])

        # gs.update(wspace=kwargs['wspace'], hspace=kwargs['hspace'],)

    else:
        fig = plt.figure(figsize=(
            (ncol+2)*kwargs["figscale"],
            nrow*kwargs["magical_factor"]*kwargs["figscale"])
        )

        gs = GridSpec(
            nrow, ncol+2,
            height_ratios=[1]*nrow,
            width_ratios=[1]*ncol + [kwargs["blankratio"], kwargs["axcbar_ratio"]])

        gs.update(wspace=kwargs['wspace'], hspace=kwargs['hspace'],)

    ims = {}
    iter_grid = [i for i in product(
        [i for i in range(nrow)],
        [j for j in range(ncol)])]

    if kwargs['offset_coordinates']:
        for n, (i, j) in enumerate(iter_grid):
            ims[n] = plt.subplot(gs[i, j])
            ax_transform = ims[n].transData
            if kwargs['box_offset'] is None:
                ims[n].set_xlim([x0y0lim[0], xNyNlim[0]])
                ims[n].set_ylim([x0y0lim[1], xNyNlim[1]])
            else:
                ims[n].set_xlim(
                    [kwargs['box_offset'][0][0],
                     kwargs['box_offset'][1][0]])
                ims[n].set_ylim(
                    [kwargs['box_offset'][0][1],
                     kwargs['box_offset'][1][1]])
            xaxis = ims[n].get_xaxis()
            yaxis = ims[n].get_yaxis()

            ims[n].tick_params(
                axis="both",
                which="both",
                direction=kwargs['tick_direction'],
                grid_color=kwargs['grid_color'],
                labelcolor=kwargs['label_color'],
                colors=kwargs['grid_color'],
                bottom=True,
                top=True,
                left=True,
                right=True,
                grid_linewidth=kwargs['grid_linewidth'],
                width=kwargs['tick_width'],
                length=kwargs['tick_length'],
                )
            ims[n].tick_params(
                axis="x",
                pad=kwargs['tick_x_pad']
            )
            ims[n].tick_params(
                axis="y",
                pad=kwargs['tick_y_pad']
            )

            if kwargs["display_x_minor_ticks"]:
                ims[n].minorticks_on()
            # xaxis.display_minor_ticks(kwargs['display_x_minor_ticks'])
            # yaxis.display_minor_ticks(kwargs['display_y_minor_ticks'])

            if kwargs["minortick_length"] == None:
                ims[n].tick_params(
                    which='minor',
                    length=kwargs['tick_length'] / 3)
                ims[n].tick_params(
                    which='minor',
                    length=kwargs['tick_length'] / 3)
            else:
                ims[n].tick_params(
                    which='minor',
                    length=kwargs['minortick_length'])
                ims[n].tick_params(
                    which='minor',
                    length=kwargs['minortick_length'])

            if kwargs['show_ticks_axis'] is None:
                if (j > 0) and (i < nrow-1):
                    xaxis.set_ticklabels([])
                    yaxis.set_ticklabels([])
                if (i == (nrow-1)) and (j > 0):
                    yaxis.set_ticklabels([])
                if (i < (nrow-1)) and (j == 0):
                    xaxis.set_ticklabels([])
            else:
                if (j == kwargs['show_ticks_axis'][1]) \
                   and (i == kwargs['show_ticks_axis'][0]):
                    pass
                else:
                    xaxis.set_ticklabels([])
                    yaxis.set_ticklabels([])

            if kwargs["tick_spacing"] is not None:
                xaxis.set_major_locator(
                    ticker.MultipleLocator(kwargs["tick_spacing"]))
                yaxis.set_major_locator(
                    ticker.MultipleLocator(kwargs["tick_spacing"]))

    else:
        for n, (i, j) in enumerate(iter_grid):
            ims[n] = plt.subplot(gs[i, j], projection=wcs)
            ims[n].set_xlim([box[0][0], box[1][0]])
            ims[n].set_ylim([box[0][1], box[1][1]])
            xaxis = ims[n].coords[0]
            yaxis = ims[n].coords[1]
            # yaxis.set_ticks([31.26638888, 31.26666667, 31.26694444,
            #                  31.26722222, 31.2675, 31.26777778] * u.deg)


            if kwargs['ticksdefaults']:
                pass
            else:
                xaxis.tick_params(direction=kwargs['tick_direction'],
                                  grid_color=kwargs['grid_color'],
                                  colors=kwargs['grid_color'],
                                  labelcolor=kwargs['label_color'],
                                  color=kwargs['grid_color'],
                                  grid_linewidth=kwargs['grid_linewidth'],
                                  width=kwargs['tick_width'],
                                  length=kwargs['tick_length'],
                                  pad=kwargs['tick_x_pad'])
                yaxis.tick_params(direction=kwargs['tick_direction'],
                                  grid_color=kwargs['grid_color'],
                                  colors=kwargs['grid_color'],
                                  labelcolor=kwargs['label_color'],
                                  color=kwargs['grid_color'],
                                  grid_linewidth=kwargs['grid_linewidth'],
                                  width=kwargs['tick_width'],
                                  length=kwargs['tick_length'],
                                  pad=kwargs['tick_y_pad'])

                xaxis.set_ticks(spacing=kwargs['tick_spacing'] * u.arcsec)
                xaxis.display_minor_ticks(kwargs['display_x_minor_ticks'])
                yaxis.display_minor_ticks(kwargs['display_y_minor_ticks'])


            if kwargs["minortick_length"] == None:
                xaxis.tick_params(which='minor',
                                  length=kwargs['tick_length'] / 3)
                yaxis.tick_params(which='minor',
                                  length=kwargs['tick_length'] / 3)

            else:
                xaxis.tick_params(which='minor', length=kwargs['minortick_length'])
                yaxis.tick_params(which='minor', length=kwargs['minortick_length'])


            # import pdb; pdb.set_trace()
            ims[n].set_xlabel(' ')
            ims[n].set_ylabel(' ')

            if kwargs['rotate_ticktext_yaxis'] is not None:
                yaxis.set_ticklabel(rotation=kwargs['rotate_ticktext_yaxis'])

            if kwargs['rotate_ticktext_xaxis'] is not None:
                xaxis.set_ticklabel(rotation=kwargs['rotate_ticktext_xaxis'])

            if kwargs['show_ticks_axis'] is None:
                if (j > 0) and (i < nrow-1):
                    xaxis.set_ticklabel_visible(False)
                    yaxis.set_ticklabel_visible(False)
                if (i == (nrow-1)) and (j > 0):
                    yaxis.set_ticklabel_visible(False)
                if (i < (nrow-1)) and (j == 0):
                    xaxis.set_ticklabel_visible(False)
            else:
                if (j == kwargs['show_ticks_axis'][1]) \
                   and (i == kwargs['show_ticks_axis'][0]):
                    pass
                else:
                    xaxis.set_ticklabel_visible(False)
                    yaxis.set_ticklabel_visible(False)

            if n in kwargs['invisible_axis']:
                ims[n].set_axis_off()

        if kwargs['dropticklabels'] is not None:
            if type(kwargs['dropticklabels']) != list:
                dtls = [kwargs['dropticklabels']]
            else:
                dtls = kwargs['dropticklabels']
            for dtl in dtls:
                if type(dtl['icoord']) != list:
                    dtl['icoord'] = [dtl['icoord']]
                else:
                    pass
                iter_drop = product(dtl['naxs_drop'], dtl['icoord'])
                for nax_drop, icoord in iter_drop:
                    ims[nax_drop].coords[icoord].ticklabels.set_dropticklabels(
                        indexes=dtl['itick'],
                        axis=dtl['axis']
                    )

    if kwargs['plot_cbar']:
        if kwargs['cbar_end'] is None:
            ax_cbar = plt.subplot(gs[:, -1])
        else:
            ax_cbar = plt.subplot(gs[:kwargs['cbar_end'], -1])

    imax = {}
    nchannels = len(image_cube)
    if kwargs['v_channels'] is None:
        v_channels = mf.vel_from_header(kwargs['header'])
    else:
        v_channels = kwargs['v_channels']
    props_dict = dict(boxstyle=kwargs['boxstyle'],
                      facecolor=kwargs['box_facecolor'],
                      alpha=kwargs['alpha_box'])
    props = props_dict if kwargs['textbox'] else None

    if kwargs['channels2plot'] is not None:
        channels2plot = kwargs['channels2plot']
    else:
        channels2plot = np.array([int(round(i)) for i in
                                  np.linspace(chan_0, chan_f, nrow*ncol)])

    for n, channel in enumerate(channels2plot):
        if kwargs['norm'] == 'linear':
            norm = colors.Normalize(vmax=kwargs['vmax'], vmin=kwargs['vmin'])
        elif kwargs['norm'] == 'log':
            norm = colors.LogNorm(vmin=kwargs['vmin'], vmax=kwargs['vmax'])
        elif kwargs['norm'] == 'symlog':
            norm = colors.SymLogNorm(linthresh=kwargs['linthresh'],
                                     linscale=kwargs['linscale'],
                                     vmin=kwargs['vmin'],
                                     vmax=kwargs['vmax'])
        elif kwargs['norm'] == 'divnorm':
            norm = colors.TwoSlopeNorm(vcenter=kwargs['vcenter'],
                                       vmin=kwargs['vmin'],
                                       vmax=kwargs['vmax'])

        image = collapse_chans(image_cube, operation, channel,
                               width_chan, box) \
            if operation is not None else image_cube[channel]
        imax[n] = ims[n].imshow(image,
                                origin=kwargs['origin'],
                                cmap=kwargs['cmap'],
                                aspect=kwargs['cube_aspect'],
                                norm=norm,
                                extent=extent,
                               # interpolation=kwargs['interpolation']
                               )
        if kwargs["ax_facecolor"] is not None:
            ims[n].set_facecolor(kwargs["ax_facecolor"])

        if kwargs['vrel']:
            vch = v_channels[channel] - kwargs['vsys']
        else:
            vch = v_channels[channel]

        vel_str = f"{vch:.{kwargs['vel_decimals']}f} {kwargs['chansunit']}"
        ims[n].text(kwargs['x_box'],
                    kwargs['y_box'],
                    vel_str,
                    fontsize=kwargs['box_fontsize'],
                    color=kwargs['textcolor'],
                    verticalalignment='top',
                    transform=ims[n].transAxes,
                    bbox=props,
                    zorder=100)

        if kwargs['add_stars']:
            if kwargs['offset_coordinates']:
                vla4b_coords = SkyCoord(*kwargs['vla4b_deg'], unit='deg')
                vla4a_coords = SkyCoord(*kwargs['vla4a_deg'], unit='deg')

                offset_frame = center_skycoords.skyoffset_frame()
                vla4a_offset = vla4a_coords.transform_to(offset_frame)

                ims[n].plot([vla4a_offset.lon.arcsec, 0],
                        [vla4a_offset.lat.arcsec, 0],
                        marker=kwargs['markerstar_style'],
                        color=kwargs['markerstar_color'],
                        linestyle='',
                        markersize=kwargs['markerstar_size'],
                        mew=kwargs['markerstar_width'],
                        zorder=kwargs['markerstar_zorder'])
            else:
                xs_star = kwargs['vla4a_deg'][0], kwargs['vla4b_deg'][0]
                ys_star = kwargs['vla4a_deg'][1], kwargs['vla4b_deg'][1]
                ims[n].plot(xs_star,
                            ys_star,
                            marker=kwargs['markerstar_style'],
                            color=kwargs['markerstar_color'],
                            linestyle='',
                            transform=ims[n].get_transform('icrs'),
                            markersize=kwargs['markerstar_size'],
                            mew=kwargs['markerstar_width'])
                if kwargs['star_nax'] is not None:
                    vla4b_coords = SkyCoord(*kwargs['vla4b_deg'], unit='deg')
                    vla4a_coords = SkyCoord(*kwargs['vla4a_deg'], unit='deg')
                    vla4b_pix = skycoord_to_pixel(vla4b_coords, wcs=wcs)
                    vla4a_pix = skycoord_to_pixel(vla4a_coords, wcs=wcs)
#                    import pdb; pdb.set_trace()
                    ims[kwargs['star_nax']].text(vla4b_pix[0]+kwargs['vla4b_offset_label_x'],
                                                 vla4b_pix[1]+kwargs['vla4b_offset_label_y'],
                                                 'B',
                                                 color='w')
                    ims[kwargs['star_nax']].text(vla4a_pix[0]+kwargs['vla4a_offset_label_x'],
                                                 vla4a_pix[1]+kwargs['vla4a_offset_label_y'],
                                                 'A',
                                                 color='w')

        if ellip_dic is not None:
            bb_types_with_chan = [bb_type for bb_type in ellip_dic if channel
                                  in ellip_dic[bb_type]['data'].index]
            beam_size_pix = (kwargs['header']['BMAJ']
                             + kwargs['header']['BMIN']) \
                * 0.5 / kwargs['header']['CDELT2']

            for bb_type in bb_types_with_chan:
                if ellip_dic[bb_type]['plot_ellipses'] == 'annulus':
                    aper = EllipticalAnnulus(
                        positions=(
                                   ellip_dic[bb_type]['data'].x0[channel],
                                   ellip_dic[bb_type]['data'].y0[channel]),
                        a_in=ellip_dic[bb_type]['data'].sma[channel]
                             - beam_size_pix / 2.,
                        a_out=ellip_dic[bb_type]['data'].sma[channel]
                              + beam_size_pix / 2.,
                        b_out=(ellip_dic[bb_type]['data'].sma[channel]
                              + beam_size_pix / 2.)
                        * (1. - ellip_dic[bb_type]['data'].eps[channel]),
                        b_in=None,
                        theta=ellip_dic[bb_type]['data'].pa[channel])

                elif ellip_dic[bb_type]['plot_ellipses'] == 'ellipse':
                    aper = EllipticalAperture((
                        ellip_dic[bb_type]['data'].x0[channel],
                        ellip_dic[bb_type]['data'].y0[channel]),
                        ellip_dic[bb_type]['data'].sma[channel],
                        ellip_dic[bb_type]['data'].sma[channel]
                        * (1. - ellip_dic[bb_type]['data'].eps[channel]),
                        ellip_dic[bb_type]['data'].pa[channel])

                elif ellip_dic[bb_type]['plot_ellipses'] == \
                        'ellipse_beam_upper':
                    aper = EllipticalAperture((
                        ellip_dic[bb_type]['data'].x0[channel],
                        ellip_dic[bb_type]['data'].y0[channel]),
                        ellip_dic[bb_type]['data'].sma[channel]
                        + beam_size_pix / 2.,
                        (ellip_dic[bb_type]['data'].sma[channel]
                            + beam_size_pix / 2.)
                        * (1. - ellip_dic[bb_type]['data'].eps[channel]),
                        ellip_dic[bb_type]['data'].pa[channel])

                aper.plot(ims[n],
                          color=ellip_dic[bb_type]['color'],
                          linewidth=ellip_dic[bb_type]['linewidth'],
                          alpha=ellip_dic[bb_type]['alpha'],
                          zorder=10000,
                          linestyle=ellip_dic[bb_type]['linestyle'])

                if channel in ellip_dic[bb_type]['chans2label']:
                    ims[n].text(ellip_dic[bb_type]['chans2label'][channel][0],
                                ellip_dic[bb_type]['chans2label'][channel][1],
                                ellip_dic[bb_type]['ring_label'],
                                color=ellip_dic[bb_type]['label_color'],
                                fontsize=ellip_dic[bb_type]['label_fontsize'])

        if kwargs['models'] is not None:
            try:
                for bb in kwargs['models']:
                    contour_levels = \
                            [np.max(kwargs['models'][bb]['data'][channel])
                             * kwargs['models'][bb]['maxfactor']]
                    show_slice_cube(
                        kwargs['models'][bb]['data'],
                        ax=ims[n],
                        channel=channel,
                        render='contours',
                        contour_colors='w',
                        contour_levels=contour_levels,
                        contour_linewidths=kwargs['contour_linewidths'],
                        plot_cbar=False,
                        plot_vel=False,
                        box=None,
                        wcs=kwargs['models'][bb]['wcs'],
                        header=kwargs['models'][bb]['hdr'],
                        add_beam=False,
                        add_scalebar=False,
                        return_ax=False,
                        font_size=kwargs['font_size'],
                        icrs_xlabel=' ',
                        icrs_ylabel=' ')

                xaxis = ims[n].coords[0]
                yaxis = ims[n].coords[1]
                xaxis.tick_params(direction=kwargs['tick_direction'],
                                  grid_color=kwargs['grid_color'],
                                  colors=kwargs['grid_color'],
                                  labelcolor=kwargs['label_color'],
                                  color=kwargs['grid_color'],
                                  grid_linewidth=kwargs['grid_linewidth'],
                                  width=kwargs['tick_width'],
                                  length=kwargs['tick_length'],
                                  pad=kwargs['tick_x_pad'])
        #        import pdb; pdb.set_trace()
                yaxis.tick_params(direction=kwargs['tick_direction'],
                                  grid_color=kwargs['grid_color'],
                                  colors=kwargs['grid_color'],
                                  labelcolor=kwargs['label_color'],
                                  color=kwargs['grid_color'],
                                  grid_linewidth=kwargs['grid_linewidth'],
                                  width=kwargs['tick_width'],
                                  length=kwargs['tick_length'],
                                  pad=kwargs['tick_y_pad'])

                xaxis.display_minor_ticks(kwargs['display_x_minor_ticks'])
                yaxis.display_minor_ticks(kwargs['display_y_minor_ticks'])

                xaxis.tick_params(which='minor', length=kwargs['tick_length'] / 3)
                yaxis.tick_params(which='minor', length=kwargs['tick_length'] / 3)

                xaxis.set_ticks(spacing=kwargs['tick_spacing'] * u.arcsec)
                ims[n].set_xlabel(' ')
                ims[n].set_ylabel(' ')

                if kwargs['rotate_ticktext_yaxis'] is not None:
                    yaxis.set_ticklabel(rotation=kwargs['rotate_ticktext_yaxis'],
                                        exclude_overlapping=True)

                if kwargs['rotate_ticktext_xaxis'] is not None:
                    xaxis.set_ticklabel(rotation=kwargs['rotate_ticktext_xaxis'],
                                        exclude_overlapping=True)

                # if you want to change the fontsize of the ticklabels, do it inside
                # show_slice_cube

            except:
                pass

    if kwargs["plot_cbar"]:
        cbar = plt.colorbar(imax[0],
                            cax=ax_cbar,
                            extend=kwargs['cbar_extend'],
                            pad=kwargs['cbar_pad'],
                            fraction=kwargs['cbar_fraction'],
                            shrink=kwargs['cbar_shrink'],
                            panchor=kwargs['cbar_panchor'],
                            anchor=kwargs['cbar_anchor'],
                            aspect=kwargs['cbar_aspect'],)
        ax_cbar.tick_params(labelsize=kwargs['font_size'])
        ax_cbar.set_yscale("linear")

        if kwargs['colorbar_units'] == 'mJy/beam':
            cbar_ticklabels = ['{:.0f}'.format(i*1000) for i in cbar.ax.get_yticks()]
            cbar.ax.set_yticklabels(cbar_ticklabels)
        else:
            pass

        # import pdb; pdb.set_trace()

        cbar.set_label(r'{}'.format(kwargs['colorbar_units']),
                       labelpad=kwargs['cbar_label_pad'],
                       fontsize=kwargs['font_size'],
                       y=kwargs['cbar_label_y'],
                       rotation=kwargs['cbar_label_rotation'],
                       )

        plt.figtext(kwargs['figtext_x_hor'], kwargs['figtext_x_vert'],
                    kwargs['icrs_xlabel'])
        plt.figtext(kwargs['figtext_y_hor'], kwargs['figtext_y_vert'],
                    kwargs['icrs_ylabel'], rotation='vertical')

    if kwargs['add_beam']:
        #  pa in radians
        pa = kwargs['header']['BPA'] * np.pi/180. + np.pi/2
        # in radians
        # semi-major axis in pixels or arcsec if offset_coordinates
        a = kwargs['header']['BMAJ'] * 3600 if kwargs['offset_coordinates'] \
                else kwargs['header']['BMAJ'] / kwargs['header']['CDELT2']
        # semi-minor axis in pixels or arcsec if offset_coordinates
        b = kwargs['header']['BMIN'] * 3600 if kwargs['offset_coordinates'] \
                else kwargs['header']['BMIN'] / kwargs['header']['CDELT2']
        if (box is None) or (kwargs['offset_coordinates']):
            xpos = kwargs['xpos_beam']
            ypos = kwargs['ypos_beam']
        else:
            xpos = box[0][0] + kwargs['xpos_beam']
            ypos = box[0][1] + kwargs['ypos_beam']
        if kwargs["offset_coordinates"]:
            geometry = EllipseGeometry(x0=xpos, y0=ypos, sma=a*0.5, eps=(1-b/a),
                                       pa=np.pi-pa)
        else:
            geometry = EllipseGeometry(x0=xpos, y0=ypos, sma=a*0.5, eps=(1-b/a),
                                       pa=pa)

        aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                                  geometry.sma*(1 - geometry.eps),
                                  geometry.pa)
        aper.plot(ims[kwargs['beam_nax']],
                  color=kwargs['beam_color'],
                  linewidth=kwargs['beam_linewidth'],
                  fill=kwargs['beam_fill'],
                  zorder=10000,
                 )
        if kwargs["add_beam_rectangle"]:
            aper_rectangle = RectangularAperture(
                (geometry.x0, geometry.y0),
                geometry.sma * 2 * kwargs['beam_factor_rectangle'],
                geometry.sma * 2 * kwargs['beam_factor_rectangle'])
            aper_rectangle.plot(ims[kwargs['beam_nax']],
                      color=kwargs['beam_color_rectangle'],
                      ec=kwargs['beam_edgecolor_rectangle'],
                      fc=kwargs['beam_color_rectangle'],
                      fill=kwargs['beam_fill_rectangle'],
                      linewidth=kwargs['beam_linewidth_rectangle'],
                      zorder=0,
                      )

    if kwargs['add_scalebar']:
        fontprops = fm.FontProperties(size=kwargs['scalebar_fontsize'])
        arcsec_per_pix = 1 if kwargs["offset_coordinates"] \
            else kwargs['header']['CDELT1'] * 3600

        scalebar_distance = (kwargs['scalebar_distance']/arcsec_per_pix) \
            * (kwargs['scalebar_units'] == 'arcsec') \
            + (kwargs['scalebar_distance']/kwargs['SVS13_distance']) \
            / arcsec_per_pix * (kwargs['scalebar_units'] == 'au')
        scalebar = AnchoredSizeBar(ims[kwargs['scalebar_nax']].transData,
                                   scalebar_distance,
                                   str(kwargs['scalebar_distance'])
                                   + ' ' + kwargs['scalebar_units'],
                                   kwargs['scalebar_loc'],
                                   pad=kwargs['scalebar_pad'],
                                   color=kwargs['scalebar_color'],
                                   frameon=False,
                                   size_vertical=kwargs['scalebar_width'],
                                   label_top=kwargs['scalebar_labeltop'],
                                   sep=kwargs['scalebar_sep'],
                                   fontproperties=fontprops)
        ims[kwargs['scalebar_nax']].add_artist(scalebar)

    if points_bluearm is not None:
        for i in range(5):
            idx_sorted = np.argsort(points_bluearm[168][:, 1])
            ims[i].plot(points_bluearm[168][idx_sorted][:, 0],
                        points_bluearm[168][idx_sorted][:, 1],
                        '--w',
                        linewidth=4,
#                        alpha=0.8,
                        transform=ims[i].get_transform('icrs'))

        for i in [i+5 for i in range(5)]:
            idx_sorted = np.argsort(points_bluearm[196][:, 1])
            ims[i].plot(points_bluearm[196][idx_sorted][:, 0],
                        points_bluearm[196][idx_sorted][:, 1],
                        '--w',
                        linewidth=4,
#                        alpha=0.5,
                        transform=ims[i].get_transform('icrs'))


    if output_name is not None:
        fig.savefig('{}{}.{}'.format(kwargs['path_save'],
                                     output_name, kwargs['output_format']),
                    bbox_inches=kwargs['bbox_inches'],
                    dpi=kwargs["dpi"]
                   )
#                    dpi=kwargs['dpi'])

    if return_ims:
        return ims


