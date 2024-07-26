from matplotlib import cm
from matplotlib import colormaps
from matplotlib import colors

import numpy as np

import os

from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy import constants as const

import subprocess

process = subprocess.Popen(["whoami"], stdout=subprocess.PIPE)
result = process.communicate()[0]
user = result.decode('utf-8').rstrip('\n')

default_params = {
    'fig_size':8,
    'origin':'lower',
    'cmap':'viridis',
    'vmin':None,
    'vmax':0.06,
    'v0':109.678,
    'vf':-93.3969,
    'vrel': False,
    'scatter_size':10,
    'scatter_alpha':0.2,
    'colorbar_shrink':1,
    'colorbar_aspect':20,
    'colorbar_ngrid':0,
    'colorbar_orientation':'vertical',
    'colorbar_fraction':0.15,
    'colorbar_pad':0.05,
    'colorbar_anchor':(0.0,0.5),
    'colorbar_panchor':(1.0,0.5),
    'colorbar_show':True,
    'colorbar_labelpad':10,
    'colorbar_ticks':None,
    'colorbar_format':'%.2f',
    'colorbar_nticks':None,
    'colorbar_invertticklabels': False,
    'colorbar_label_SMA_dist': 'Deprojected Velocity (km/s)',
    'colorbar_negative_ticklabels': True,
    'colorbar_outlinelinewidth': 0.5,
    'colorbar_units': 'Jy/beam',
    'colorbar_strformat': None,
    'plotregr': True,
    'cbar_end': None,
    'figscale': 2,
    'blankratio': 0.1,
    'axcbar_ratio': 0.1,
    'cbar_label_y': 0.5,
    'cbar_label_rotation': 90,
    'invisible_axis': [],
    'x_label':'x',
    'y_label':'y',
    'z_label':'z',
    'textbox':None,
    'textcolor':'w',
    'boxstyle':'round',
    'box_facecolor':'white',
    'facecolor_background':'w',
    'plot_vel':True,
    'vel_decimals': 1,
    'v_channels': None,
    'cbar_unit':'Jy/beam',
    'offset_coordinates': False,
    'render':'raster',
    'contour_colors':'k',
    'contour_levels':5,
    'contour_linewidths':1,
    'contour_sigma_filter':None,
    'contour_zorder': 2,
    'contour_alpha': 1,
    'overplot_contours': True,
    'plot_cbar':True,
    'x_box': 0.05,
    'y_box': 0.95,
    'box_fontsize': 15,
    'box_offset': None,
    'fraction':0.1,
    'pad':0.1,
    'save_fig':False,
    'show_plot':True,
    'save_format':'eps',
    'interpolation':None,
    'filterrad':4,
    'norm':'linear',
    'vcenter':0.05,
    'linthresh':0.03,
    'linscale':0.03,
    'wcs':None,
    'cube_aspect':'equal',
    'add_stars':True,
    'star_nax': None,
    'vla4b_offset_label_x': -7,
    'vla4b_offset_label_y': 15,
    'vla4a_offset_label_x': -7,
    'vla4a_offset_label_y': 15,
    'markerstar_size':10,
    'markerstar_zorder': 5,
    'markerstar_style':'+',
    'markerstar_color':'w',
    'markerstar_width':1,
    'markerbluearm_style':'.',
    'cbar_label_pad':20,
    'cbarvmin': None,
    'cbarvmax': None,
    'label_cbar': "",
    'cbar_pad':0.05,
    'cbar_fraction':0.1,
    'cbar_extend':'both',
    'cbar_shrink':1,
    'cbar_anchor':(0.0, 0.5),
    'cbar_panchor':(1.0, 0.5),
    'cbar_aspect':20,
    'cbar_setticks': None,
    'icrs_xlabel':'ICRS RA',
    'icrs_ylabel':'ICRS DEC',
    'xlabelpad': 1,
    'ylabelpad': -0.5,
    'pv_xlabel':'Offset (arcsec)',
    'pv_ylabel': r'$V_{\mathrm{LSR}}$(km/s)',
    'rotate_ticktext_yaxis': None,
    'rotate_ticktext_xaxis': None,
    'display_x_minor_ticks': False,
    'display_y_minor_ticks': False,
    'dropticklabels': None,
    'alpha_imshow': None,
    'font_size':20,
    'alpha_box':0.85,
    'magical_factor':9.7,
    'wspace':0.01,
    'hspace':0.01,
    'tick_spacing':1.5, #arcsec
    'tick_y_spacing':1, #arcsec
    'tick_x_spacing':1.5, #arcsec
    'tick_direction':'in',
    'grid_color':'w',
    'label_color':'k',
    'subgridspec': None,
    'grid_linewidth':10,
    'tick_width':1,
    'tick_length':10,
    'minortick_length': 3,
    'tick_x_pad': 1,
    'tick_y_pad': 1,
    'show_ticks_axis': None,
    'output_format':'eps',
    'bbox_inches':'tight',
    'path_save':'/home/{}/radio_astronomy/SVS13/paperplots/'.format(user),
    'path_save_bb_characterization':'/home/{}/radio_astronomy/SVS13/bb_characterization/'.format(user),
    'path_video':'/home/{}/radio_astronomy/SVS13/paperplots/video/'.format(user),
    'path_bb_data':'/home/{}/radio_astronomy/SVS13/burbujas_data/'.format(user),
    'path_database':'/home/{}/radio_astronomy/SVS13/databases/'.format(user),
#    'path_fits':'/home/{}/data/SVS13_nocont.fits'.format(user),
#     'path_fits':'/home/{}/data/spw-2-9-gb.contsub.allchans.subimage.fits'.format(user),
#    'path_fits':'/home/{}/data/spw-2-9-gb.contsub.hanning.lsr2.fits'.format(user),
    'path_folder_points':'/home/{}/radio_astronomy/SVS13/regions_arms/spw-2-9-gb.contsub.lsr.'.format(user),
    'path_streamer_points':
        '/home/{}/radio_astronomy/SVS13/streamers/regions_streamers'.format(user),
    'path_bs': '/home/{}/radio_astronomy/SVS13/bowshock'.format(user),
    'dpi': 150,
    'use_tex':False,
    'header':None,
    'add_beam':False,
    'beam_linewidth':1.5,
    'beam_color':'w',
    'beam_fill': False,
    'beam_factor_rectangle': 1.5,
    'add_beam_rectangle': False,
    'beam_color_rectangle':'w',
    'beam_edgecolor_rectangle':'k',
    'beam_fill_rectangle': False,
    'beam_nax':0,
    'xpos_beam':30,
    'ypos_beam':30,
    'add_scalebar':False,
    'scalebar_distance':100,
    'scalebar_fontsize':15,
    'scalebar_width':1.,
    'scalebar_color':'white',
    'scalebar_loc':'lower right',
    'scalebar_pad':1,
    'scalebar_units':'au',
    'scalebar_nax':0,
    'scalebar_labeltop':False,
    'scalebar_sep':5.,
    'scalebar_bbox_to_anchor': None,
    'add_scalebar_rectangle': False,
    "scalebar_coords_rectangle": (0,0),
    "scalebar_width_rectangle": 1,
    "scalebar_height_rectangle": 1,
    'scalebar_color_rectangle': "k",
    'scalebar_edgecolor_rectangle': "k",
    'scalebar_color_rectangle': "k",
    'scalebar_fill_rectangle': True,
    'SVS13_distance':300, #pc
    'SVS13_vsys':8.5,
    'vla4b_vsys': 9.33,
    'figtext_x_hor':0.4,
    'figtext_x_vert':0.08,
    'figtext_y_hor':0.04,
    'figtext_y_vert':0.5,
    'n_sigma':3,
    'cmap_bg':'gray',
    'cmap_contours':'coolwarm_r',
    'cmap_nanized':'cool_r',
    'cmap_ellipses':'viridis',
    'ellipse_linewidth': 1,
    'plot_contours':True,
    'contour_area':True,
    'contour_alpha':0.25,
    'contour_levels':0,
    'contour_styles': 'solid',
    'theta1_pa':0,
    'theta2_pa':360,
    'show_slice_return':None,
    'bb_centers':None,
    'markercenter_color':'',
    'markercenter_style':'*',
    'markercenter_size':2,
    'markercenter_width':0.1,
    'ellipse_color':'r',
    'markerellfit_color':'b',
    'markerellfit_style':'.',
    'markerellfit_size':0.5,
    'hodapp_zorder':6,
    'stars_zorder':7,
    'color_ellip_results':'r',
    'linewidth_ellip_results':2,
    'av_addjetdir':True,
    'av_linecolor':'r',
    'av_linewidth':2,
    'av_linestyle':'solid',
    'pv_vmax':0.08,
    'vcenter_factor':1/3.,
    'pvline_color':'b',
    'pvline_width':1,
    'pvline_style':'--',
    'pvline_trans_color':'r',
    'pvline_trans_width':1,
    'pvline_trans_style':'--',
    'pvtext_x':-1.25,
    'pvtext_y':-85,
    'pv_jetline_x':None,
    'pv_jetlinecolor':'r',
    'pv_jetlinewidth':2,
    'pv_jetlinestyle':'solid',
    'twosided_plot':False,
    'channels2plot':None,
    'chansunit': "km/s",
    'inverse_plot':False,
    'plot_legend':False,
    'n_points_ellipse':50,
    'n_points_PV':50,
    'text_x_pos': [600, 730],
    'text_y_pos': [860, 1080],
    'plot_special_points': True,
    'plot_transversal_pv': True,
    'plot_ellipse_fit': True,
    'models': None,
}

def list2str(a, precision=2):
    _list = [float(f'{i:.{precision}f}') for i in a]
    _str = str(_list) if len(_list)>1 else str(_list[0])
    return _str

def get_color(vel_range, vel, cmap, norm="linear"):
    """
    Gets the color that corresponds in a colormap linearly interpolated taking
    into account the values at the limits.
    """
    cmapp = colormaps.get_cmap(cmap)
    if norm == "linear":
        norm = colors.Normalize(vmin=vel_range[0], vmax=vel_range[-1])
    elif norm == "log":
        norm = colors.LogNorm(vmin=vel_range[0], vmax=vel_range[-1])
    rgba = cmapp(norm(vel))
    color = colors.to_hex(rgba)
    return color

def make_sequential_cm(ref_cm_name, n_seq=20, n_ref=20):
    ref_cm = cm.get_cmap(ref_cm_name, n_ref)(np.linspace(0, 1, n_ref))
    seq_cmap = {}
    for i, color in enumerate(ref_cm):
        seq_cmap[i] = [list(color[:3]) + [alpha]
                       for alpha in list(np.linspace(0, 1, n_seq))]
    return seq_cmap


def get_ak_bbcenters():
    nbb = [i+1 for i in range(4)]
    bb_pos_path = {i: default_params['path_bb_data']+'pos_deg_{}.dat'.format(i)
                   for i in nbb}
    bb_files = {i: open(bb_pos_path[i]) for i in nbb}
    make_list_float = lambda l: [float(x) for x in l]
    bb_pos = {i: np.array([make_list_float(line.split('\t')[1:3])
                           for k, line in enumerate(bb_files[i]) if k > 0])
              for i in nbb}
    return bb_pos


def vel_from_header(header, vel_axis=''):
    # if header['CTYPE3'] == 'FREQ':
    #     rest_freq = header['RESTFRQ'] * u.Hz
    #     hz_per_chan = header['CDELT3'] * u.Hz
    #     chan0_freq = header['CRVAL3'] * u.Hz
    #     vel_per_chan = (const.c / rest_freq * hz_per_chan).to(u.km/u.s)
    #     v0 = (const.c / rest_freq * (rest_freq-chan0_freq)).to(u.km/u.s)
    #     vf = (v0 - vel_per_chan * (header['NAXIS3']-1)).to(u.km/u.s)
    #     v_channels = np.linspace(v0.value, vf.value, header['NAXIS3'])
    # 2021-02-23: Some modifications done for reading cut images with
    # imsubimage casa task:
    if header['CTYPE3'] == 'FREQ':
        rest_freq = header['RESTFRQ'] * u.Hz
        hz_per_chan = header['CDELT3'] * u.Hz
        chan0_freq = header['CRVAL3'] * u.Hz
        vel_per_chan = (const.c / rest_freq * hz_per_chan).to(u.km/u.s)
        v0 = (const.c / rest_freq * (rest_freq-chan0_freq)).to(u.km/u.s)
        vf = (v0 - vel_per_chan *
              (header['NAXIS3']-header["CRPIX3"])).to(u.km/u.s)
        v_channels = np.linspace(
            v0.value,
            vf.value,
            int(header['NAXIS3']-header["CRPIX3"]+1))[-header["NAXIS3"]:]

    elif header['CTYPE3'] == 'VRAD':
        nchannels = header["NAXIS3"]
        vel_per_chan = header['CDELT3']
        v0 = header['CRVAL3']
        vf = v0 + vel_per_chan * (nchannels-1)
        v_channels = np.linspace(v0, vf, nchannels)
    elif header['CTYPE3'] == 'Dynamical Time':
        t0 = header['CRVAL3']
        tf = (t0 + header['CDELT3'] * (header['NAXIS3']-1))
        v_channels = np.linspace(t0, tf, header['NAXIS3'])

    elif header['CTYPE2'] == 'FREQ':
        rest_freq = header['RESTFRQ'] * u.Hz
        hz_per_chan = header['CDELT2'] * u.Hz
        chan0_freq = header['CRVAL2'] * u.Hz
        vel_per_chan = (const.c / rest_freq * hz_per_chan).to(u.km/u.s)
        v0 = (const.c / rest_freq * (rest_freq-chan0_freq)).to(u.km/u.s)
        vf = (v0 - vel_per_chan * (header['NAXIS2']-1)).to(u.km/u.s)
        v_channels = np.linspace(v0.value, vf.value, header['NAXIS2'])

    else:
        rest_freq = header['RESTFRQ'] * u.Hz
        hz_per_chan = header['CDELT{}'.format(vel_axis)] * u.Hz
        chan0_freq = header['CRVAL{}'.format(vel_axis)] * u.Hz
        vel_per_chan = (const.c / rest_freq * hz_per_chan).to(u.km/u.s)
        v0 = (const.c / rest_freq * (rest_freq-chan0_freq)).to(u.km/u.s)
        vf = (v0 - vel_per_chan *
              (header['NAXIS{}'.format(vel_axis)]-1)).to(u.km/u.s)
        v_channels = np.linspace(v0.value,
                                 vf.value,
                                 header['NAXIS{}'.format(vel_axis)])
    return v_channels

def freq_from_header(header):
    if header['CTYPE3'] == 'FREQ':
        # rest_freq = header['RESTFRQ'] * u.Hz
        hz_per_chan = header['CDELT3'] * u.Hz
        chan0_freq = header['CRVAL3'] * u.Hz
        refpix = header["CRPIX3"]
        nchans = header["NAXIS3"]
        f0 = chan0_freq - hz_per_chan * (refpix-1)
        ff = chan0_freq + hz_per_chan * (nchans-refpix)

        f_channels = np.linspace(
            f0.value,
            ff.value,
            int(nchans-refpix+1))[-nchans:]
        return f_channels

def chan_from_vel(vel, header):
    vs = [abs(v-vel) for v in vel_from_header(header)]
    chan_vel = np.where(vs == np.min(vs))
    return chan_vel

def change_velrange_header(header_new, header_old, chan_0, vel_axis='2'):
    freq_offset = header_old['CDELT{}'.format(vel_axis)] * chan_0
    header_new['CRVAL{}'.format(vel_axis)] = \
        header_old['CRVAL{}'.format(vel_axis)] + freq_offset


hms_to_deg = lambda h, m, s: 360*(h/24+m/(24*60)+s/(24*3600) )


def deg_to_hms(deg):
    h = int(deg*(24/360))
    h_rest = deg*(24/360)-h
    m = int(h_rest*60)
    m_rest = h_rest*60-m
    s = m_rest*60
    return h, m, s

# SkyCoord(RA*u.deg, DEC*u.deg, frame='icrs')
# c = SkyCoord(ra='03h29m3.7454s', dec='31d16m3.784s', frame='icrs')
# p.ra.hms
# p.dec


def arcs2skycoord(arcs_pix, header):
    wcs = WCS(header).celestial
    arcs_sky = {arc: {} for arc in arcs_pix}
    for arc in arcs_pix:
        arcs_sky[arc]['x0'] = pixel_to_skycoord(arcs_pix[arc]['x0'],
                                                arcs_pix[arc]['y0'],
                                                wcs).ra
        arcs_sky[arc]['y0'] = pixel_to_skycoord(arcs_pix[arc]['x0'],
                                                arcs_pix[arc]['y0'],
                                                wcs).dec
        arcs_sky[arc]['width'] = arcs_pix[arc]['width'] \
            * header['CDELT2'] * 3600 * u.arcsec
        arcs_sky[arc]['height'] = arcs_pix[arc]['height'] \
            * header['CDELT2'] * 3600 * u.arcsec
        arcs_sky[arc]['angle'] = arcs_pix[arc]['angle']
        arcs_sky[arc]['theta1'] = arcs_pix[arc]['theta1']
        arcs_sky[arc]['theta2'] = arcs_pix[arc]['theta2']
    return arcs_sky


def arcs2pix(arcs_sky, header, key_cdelt=["CDELT1", "CDELT2"]):
    wcs = WCS(header).celestial
    arcs_pix = {arc: {} for arc in arcs_sky}

    for arc in arcs_sky:
        arcs_pix[arc]['x0'] = skycoord_to_pixel(SkyCoord(arcs_sky[arc]['x0'],
                                                         arcs_sky[arc]['y0']),
                                                wcs)[0]
        arcs_pix[arc]['y0'] = skycoord_to_pixel(SkyCoord(arcs_sky[arc]['x0'],
                                                         arcs_sky[arc]['y0']),
                                                wcs)[1]
        arcs_pix[arc]['width'] = np.abs((arcs_sky[arc]['width'].to(u.deg) /
                                  (header[key_cdelt[0]] * u.deg)).value)
        arcs_pix[arc]['height'] = np.abs((arcs_sky[arc]['height'].to(u.deg) /
                                   (header[key_cdelt[1]] * u.deg)).value)
        params_arc = ['angle', 'theta1', 'theta2']

        for param in params_arc:
            arcs_pix[arc][param] = arcs_sky[arc][param] \
                if param in arcs_sky[arc] else None
            arcs_pix[arc][param] = arcs_sky[arc][param] \
                if param in arcs_sky[arc] else None
            arcs_pix[arc][param] = arcs_sky[arc][param] \
                if param in arcs_sky[arc] else None
    return arcs_pix

def arcs2offset(arcs_sky, center_coords):
    arcs_offset = {arc: {} for arc in arcs_sky}
    center_frame = center_coords.skyoffset_frame()

    for arc in arcs_sky:
        coord_sky = SkyCoord(arcs_sky[arc]['x0'],
                             arcs_sky[arc]['y0'], unit="deg")
        coord_offset = coord_sky.transform_to(center_frame)
        arcs_offset[arc]['x0'] = coord_offset.lon.to(u.arcsec).value
        arcs_offset[arc]['y0'] = coord_offset.lat.to(u.arcsec).value
        arcs_offset[arc]['width'] = arcs_sky[arc]['width'].to(u.arcsec).value
        arcs_offset[arc]['height'] = arcs_sky[arc]['height'].to(u.arcsec).value
        if "angle" in arcs_sky[arc]:
            arcs_offset[arc]["angle"] = - arcs_sky[arc]["angle"]
        if "theta1" in arcs_sky[arc]:
            arcs_offset[arc]["theta1"] = 180 - arcs_sky[arc]["theta2"]
        if "theta2" in arcs_sky[arc]:
            arcs_offset[arc]["theta2"] = 180 - arcs_sky[arc]["theta1"]
    return arcs_offset


def offset_to_sky(offset_coord, reference_skycoord, unit=u.arcsec):
    """
    Computes the skycoordinates given the offset position in arcsec of the
    target and the skycoordinates (SkyCoord type)
    """
    sep = np.sqrt(offset_coord[0]**2 + offset_coord[1]**2)
    pa = np.arctan(offset_coord[1] / offset_coord[0]) * 180 / np.pi + 90
    skyco = reference_skycoord.directional_offset_by(pa*u.deg, sep*unit)
    return skyco


def hdr_for_aper(hdr2copy, NAXIS1, NAXIS2, CRPIX1, CRPIX2, CRVAL1, CRVAL2):
    """
    Creates the header for a subimage of hdr2copy, recentering as desired
    """
    hdr_aper = hdr2copy.copy()
    remove_from_header = ["NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3",
                          "NAXIS4", "CTYPE4", "CRVAL4", "CDELT4", "CRPIX4", "CUNIT4",
                          "PC1_3", "PC1_4", "PC2_3", "PC2_4",
                          "PC3_1", "PC3_2", "PC3_3", "PC3_4",
                          "PC4_1", "PC4_2", "PC4_3", "PC4_4"]
    for key in remove_from_header:
        try:
            hdr_aper.remove(key)
        except:
            pass

    hdr_aper["NAXIS"] = 2
    hdr_aper["NAXIS1"] = NAXIS1
    hdr_aper["NAXIS2"] = NAXIS2
    hdr_aper["CRPIX1"] = CRPIX1
    hdr_aper["CRPIX2"] = CRPIX2
    hdr_aper["CRVAL1"] = CRVAL1
    hdr_aper["CRVAL2"] = CRVAL2
    return hdr_aper

def wcs_for_aper(image_original, wcs_original, center_pix, size_pix):
    """
    Gets the wcs for an aperture cutting the original image.
    """
    _cutout = Cutout2D(image_original,
                          position=tuple(center_pix),
                          size=np.array(size_pix),
                          wcs=wcs_original)
    return _cutout.wcs

def negated_array(data):
    """
    Returns the complementary of an array of zeros and ones.
    """
    array = np.copy(data)
    array_0s = array == 0
    array_1s = array == 1
    array[array_0s] = 1
    array[array_1s] = 0
    return array

def padmask(data, mask):
    """
    Given the 2D array data and a photutils.aperture.mask.ApertureMask, it
    returns a 2D array of 0's and 1's, but with the dimensions of the initial
    data, with 0 outside the mask.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The 2D data array to which the mask is applied.

    mask : `photutils.aperture.mask.ApertureMask`
        The mask region

    Returns
    ----------
    padded : 2D `~numpy.ndarray`
        The 2D data, conserving its dimensions pixel coordinates, with 1's for pixel
        within the mask and 0's for pixels outside the mask
    """
    datany, datanx = np.shape(data)
    xmin = mask.bbox.ixmin
    xmax = mask.bbox.ixmax
    ymin = mask.bbox.iymin
    ymax = mask.bbox.iymax
    mask_changed = np.copy(mask)
    if datanx >= xmax:
        padx_right = datanx - xmax
    else:
        mask_changed = mask_changed[:, :-int(np.abs(datanx-xmax))]
        padx_left = xmin #+ np.abs(datanx-xmax)
        padx_right = 0
    if datany >= ymax:
        pady_upper = datany - ymax
    else:
        mask_changed = mask_changed[:-int(np.abs(datany-ymax)), :]
        pady_lower = ymin #+ np.abs(datany-ymax)
        pady_upper = 0
    if xmin < 0:
        mask_changed = mask_changed[:, np.abs(xmin):]
        padx_left = 0
    else:
        padx_left = xmin
    if ymin < 0:
        mask_changed = mask_changed[np.abs(ymin):, :]
        pady_lower = 0
    else:
        pady_lower = ymin
    padded = np.pad(mask_changed, ((pady_lower, pady_upper), (padx_left, padx_right)))
    return padded

def cutout3D(data, chans, position, size, header, cutout2d_params={}):
    """
    Create a cutout object from a 3D array, representing a spectral cube.

    This function uses Cutout2D function from astropy.nddata

    If a `~astropy.wcs.WCS` object is input, then the returned object
    will also contain a copy of the original WCS, but updated for the
    cutout array.

    .. warning::

        The cutout WCS object does not currently handle cases where the
        input WCS object contains distortion lookup tables described in
        the `FITS WCS distortion paper
        <https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf>`__.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The 3D data array from which to extract the cutout array.

    chans: `~numpy.ndarray`
        Subset of spectral channels to cutout

    position : tuple or list, `~astropy.coordinates.SkyCoord`
        The position of the cutout array's center with respect to
        the ``data`` array.  The position can be specified either as
        a ``(x, y)`` tuple of pixel coordinates or a
        `~astropy.coordinates.SkyCoord`, in which case ``wcs`` is a
        required input.

    size : int, array_like, or `~astropy.units.Quantity`
        The size of the cutout array along each axis.  If ``size``
        is a scalar number or a scalar `~astropy.units.Quantity`,
        then a square cutout of ``size`` will be created.  If
        ``size`` has two elements, they should be in ``(ny, nx)``
        order.  Scalar numbers in ``size`` are assumed to be in
        units of pixels.  ``size`` can also be a
        `~astropy.units.Quantity` object or contain
        `~astropy.units.Quantity` objects.  Such
        `~astropy.units.Quantity` objects must be in pixel or
        angular units.  For all cases, ``size`` will be converted to
        an integer number of pixels, rounding the the nearest
        integer.  See the ``mode`` keyword for additional details on
        the final cutout size.

        .. note::
            If ``size`` is in angular units, the cutout size is
            converted to pixels using the pixel scales along each
            axis of the image at the ``CRPIX`` location.  Projection
            and other non-linear distortions are not taken into
            account.

    header : `~astropy.io.fits.header.Header`
        Header of the original data

    cutout2d_params: dict, optional
        Dictionary containing the optional parameters to pass to Cutout2D

    Returns
    ----------
    cutoutdata : 3D `~numpy.ndarray`
        The 3D cutout array.

    cutoutheader : `~astropy.io.fits.header.Header`
        A header object associated with the cutout array
    """
    chans2cut = np.arange(chans[0], chans[1]+1)
    cutout_data = np.zeros([len(chans2cut), size[0], size[1]])
    wcs = WCS(header).celestial
    for n, chan in enumerate(chans2cut):
        _cutout = Cutout2D(data[chan],
                              position,
                              size,
                              wcs=wcs,
                              **cutout2d_params)
        cutout_data[n] = _cutout.data

    hdr_cutout_2D = _cutout.wcs.to_header()
    hdr_cutout = header.copy()
    remove_from_header = [
         "NAXIS4", "CTYPE4", "CRVAL4",
         "CDELT4", "CRPIX4", "CUNIT4",
         "PC1_4", "PC2_4", "PC4_1",
         "PC4_2", "PC4_3", "PC4_4", "PC3_4",
    ]

    for key in remove_from_header:
        try:
            hdr_cutout.remove(key)
        except:
            pass

    hdr_cutout["NAXIS"] = 3
    hdr_cutout["NAXIS1"] = size[0]
    hdr_cutout["NAXIS2"] = size[1]
    hdr_cutout["NAXIS3"] = len(chans2cut)
    hdr_cutout["CRPIX1"] = hdr_cutout_2D["CRPIX1"]
    hdr_cutout["CRPIX2"] = hdr_cutout_2D["CRPIX2"]
    hdr_cutout["CRPIX3"] = 1.000000000000E+00
    hdr_cutout["CRVAL1"] = hdr_cutout_2D["CRVAL1"]
    hdr_cutout["CRVAL2"] = hdr_cutout_2D["CRVAL2"]
    hdr_cutout["CRVAL3"] = header["CRVAL3"] + header["CDELT3"]*chans[0]

    return cutout_data, hdr_cutout


def merge3D(cubes, headers, intersop="sum", allow_diffspecbinning=False):
    """
    Merges spectral cubes into one. The cubes must have the same dimensions in
    pixels and sky coordinates; the merge is only performed along the spectral
    axis, and the channel width should be the same.

    Parameters
    ----------
    cubes : list
        List of the cubes

    headers : list
        List of the headers

    intersop : str
        Operation to performed for the pixel intersection

    allow_diffspecbinning : bool
        If the frequency binning is not the same, tries to match the bins
        to the nearest. This would work ok if the binning is shifted rather
        than if the channel width is different.

    Returns
    ----------
    merged_cube : 3D `~numpy.ndarray`
        The 3D merged data.

    merged_header : `~astropy.io.fits.header.Header`
        A header associated to the 3D cube
    """

    keys = list(cubes)
    keyref = keys[0]

    freq_cubes = {key: freq_from_header(headers[key])
                            for key in keys}
    freq_set = set(np.concatenate([freq_cubes[key] for key in keys]).ravel())
    minfreq = np.min(list(freq_set))
    maxfreq = np.max(list(freq_set))
    sizey, sizex = np.shape(cubes[keyref])[1:3]

    freq_merged = np.arange(minfreq,
                            maxfreq+headers[keyref]["CDELT3"],
                            headers[keyref]["CDELT3"])

    nchans_merged = len(freq_merged)
    merged_header = headers[keyref].copy()
    merged_header["NAXIS3"] = nchans_merged
    merged_header["CRPIX3"] = 1.000000000000E+00
    merged_header["CRVAL3"] = minfreq

    chans_merged = np.arange(0, nchans_merged)

    cubes_merge_dim = {key: np.zeros((nchans_merged, sizey, sizex)) for key in keys}
    for key in keys:
        for chanmerge, freq in enumerate(freq_merged):
            if allow_diffspecbinning:
                min_fdiff = np.min(np.abs(freq_cubes[key]-freq))
                freqwidth = np.abs(freq_cubes[key][1]-freq_cubes[key][0])
                if freqwidth/2. >= min_fdiff:
                    chan = np.argmin(np.abs(freq_cubes[key]-freq))
                    # import pdb; pdb.set_trace()
                    cubes_merge_dim[key][chanmerge] = cubes[key][chan]
            else:
                if freq in freq_cubes[key]:
                    chan = np.where(freq_cubes[key]==freq)[0][0]
                    cubes_merge_dim[key][chanmerge] = cubes[key][chan]

    if intersop == "sum":
        merged_cube = np.sum([cubes_merge_dim[key] for key in keys], axis=0)
    else:
        # TODO: Implement other operations
        pass

    return merged_cube, merged_header


def get_nmax(data, n):
    data_nomax = np.copy(data)
    maxes = []
    for i in range(n):
        maxes.append(np.max(data_nomax))
        data_nomax[np.argmax(data_nomax)] = 0
    return maxes

def klambda2MRS(kl, ff):
    MRS = ff * 10**(-3) / kl * u.radian
    return MRS.to(u.arcsec)

def MRS2klambda(MRS, ff):
    kl = ff * 10**(-3) / MRS.to(u.radian).value
    return kl

def klambda2bl(kl, wl):
    bl = 10**3 * wl * kl
    return bl.to(u.m)

def bl2klambda(bl, wl):
    kl = bl.to(u.m).value / (wl.to(u.mm).value)
    return kl

def dv_func(nu0, dnu, radio_def=True):
    """
    radio_def True for radio definition, False for optical definition, which
    are both approximations:
    https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php
    Note that radio defition is deprecated by the IAU.
    """
    if radio_def:
        dv = (const.c / nu0 * (dnu)).to(u.km/u.s)
    else:
        nu = nu0 - dnu
        dv = (const.c * nu0 * (1 / nu - 1 / nu0)).to(u.km/u.s)
        # According to astropy should be:
        # dv = const.c * dnu / nu
    return dv

def dnu_func(nu0, dv, radio_def=True):
    """
    radio_def True for radio definition, False for optical definition, which
    are both approximations:
    https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php
    Note that radio definition is deprecated by the IAU.
    """
    if radio_def:
        nu = (dv * nu0 / const.c).to(u.kHz)
    else:
        nu = (nu0 / (dv / const.c + 1)).to(u.kHz)
    return nu

def doppler_shift(nu0, nu):
    """
    Calculates the doppler velocity shift
    """
    v = const.c * (nu0**2 - nu**2) / (nu0**2 + nu**2)
    return v.to(u.km/u.s)

def primary_beam_func(l, diameter):
    """
    Calculates the primary beam size approximating it to a Gaussian FHWM

    Parameters
    ----------
    l: astropy.units
        wavelength
    diamter: astropy.units
        diameter of a single dish
    """
    return (1.13 * (l/diameter) * 180 / np.pi * 3600 * u.arcsec).to(u.arcsec)

def proper_motion_func(v_los, d, i, v_sys):
    """
    Calculates the proper motion in arcsec per year of an object

    Parameters
    ---------
    v_los: astropy.units
        line of sight velocity in km/s
    d: astropy.units
        distance to the source
    i: float
        inclination angle in degrees
    v_sys: astropy.units
        systemic velocity of the source
    """
    v_ps = (v_los - v_sys) * np.tan(i * np.pi / 180)
    v_ps_arcs_yr = (v_ps.to(u.au / u.s) / d.to(u.pc) ).value * 3600 * 24 * 365
    return v_ps_arcs_yr * u.arcsec / u.yr

def Bnu_f(nu, T):
    """
    Planck's law. Intensity dimensions (erg s-1 cm-2 Hz-1 sr-1, or Jy sr-1)
    """
    Bnu = 2 * const.h * nu**3 / const.c**2 * (np.exp(const.h*nu/(const.k_B*T))-1)**(-1)
    return Bnu.to(u.Jy) / u.sr

def Tbb(nu, Inu):
    """
    Temperature of a blackbody with intensity (Jy/sr) of Inu
    """
    T = (const.h * nu / const.k_B) \
        / np.log(2*const.h*nu**3/(const.c**2*Inu*u.sr) +1)
    return T

def Inu_f(TR, nu):
    """
    Intensity as a function of radiation temperature TR
    """
    Inu = 2 * const.k_B * nu**2 * TR / const.c**2
    return Inu.to(u.Jy) / u.sr

def Jnu_f(nu, T):
    """
    Intensity in temperature units. In the R-J limit, this is Tb,
    the brightness temperature
    """
    Jnu = const.h * nu / const.k_B \
          * (np.exp(const.h*nu / (const.k_B*T)) - 1)**(-1)
    return Jnu.to(u.K)

def TR_f(nu, Tb):
    """
    Radiation temperature as a function of brightness temperature.
    For the R-J limit, TR=Tb
    """
    TR = Jnu_f(nu, Tb)
    return TR.to(u.K)

def mb_sa_gaussian_f(maja, mina):
    """
    Solid angle of a gaussian main beam and θmaj and θmin as
    the half-power beam widths
    """
    omega_M = np.pi * maja * mina / (4 * np.log(2))
    return omega_M.to(u.sr)

def Snu_obs_f(nu, TR_obs, beam_sa):
    """
    Flux density given de radiation temperature (brightness temperature
    in RJ limit) observed, and the solid angle of the main beam
    """
    Snu_obs = 2 * const.k_B * nu**2 * TR_obs * beam_sa.to(u.sr).value / const.c**2
    return Snu_obs.to(u.Jy)

def Snu_f(nu, Tex):
    """
    Flux density as a function of the excitation temperature
    """
    Snu = 2 * const.k_B * nu**2 * Jnu_f(nu, Tex) / const.c**2
    return Snu.to(u.Jy)

def T0_f(nu, Tex, tau_0, Tbg=2.7*u.K):
    """
    Transport eq. for line intensity, maximum intensity or center intesity (in temperature
    units). Tau_0 is the maximum optical depth, or optical depth in the center of the line
    """
    T0 = (Jnu_f(nu, Tex) - Jnu_f(nu, Tbg)) * (1 - np.exp(-tau_0))
    return T0.to(u.K)

def Tex_f(nu, T0, tau_0, Tbg=2.7*u.K):
    """
    Excitation temperature according to transport eq., as a function of tau_0 and
    the observed T0 (intensity in temperature units)
    """
    Tex = const.h * nu / const.k_B \
          * np.log((const.h*nu/const.k_B)*(T0/(1-np.exp(-tau_0)) + Jnu_f(nu,Tbg))**(-1) + 1)**(-1)
    return Tex.to(u.K)

def tau_0_f(nu, T0, Tex, Tbg=2.7*u.K):
    """
    Optical depth in the center of the line as a function of the intensity in temperature
    units (T0)
    """
    tau_0 = - np.log(1 - T0 / (Jnu_f(nu,Tex) - Jnu_f(nu,Tbg)))
    return tau_0

def Snu_obs_2_Inu(Snu_per_beam, beam_sa):
    """
    Calculates the intensity given the observed density flux per beam
    """
    snu = Snu_per_beam / beam_sa
    return snu.to(u.Jy / u.sr)

def Inu_2_T0(nu, Inu):
    """
    Calculates the intensity in temperature units given the intensity
    """
    T0 = (Inu*u.sr) * const.c**2 / (2 * const.k_B * nu**2)
    return T0.to(u.K)

def nu_obs_f(v_r, nu_0):
    """
    Calculates the observed frequency of a doppler shifted emission
    """
    nu_obs = nu_0 * (1 - v_r / const.c)
    return nu_obs.to(u.GHz)

def A_j_jm1_f(J, nu, dip_mom):
    """
    Calculates the spontaneus emission coeffitient for the J -> J-1 transition
    of a molecule with a dipole moment mu.
    """
    aj_jm1 = 1.165 * 10**(-11) * dip_mom.to(u.D).value**2 * (J/(2*J+1)) * nu.to(u.GHz).value**3
    return aj_jm1 * u.s**(-1)

#    aj_jm1 = 64 * np.pi**4 / (3*const.h*const.c**3) * (dip_mom**2*J/(2*J+1)) * nu**3
#    return aj_jm1

def tau_0_cd_f(A, J, nu, dv, Tex, N):
    """
    Optical depth as a function of column density N for a CO
    rotational transition  J->J-1
    Parameters:
    A: Einstein emission coeffitient
    J: Transition J->J-1
    nu: Frequency of the transition
    dv: line FWHM in velocity units
    Tex: Excitation temperature
    N: Column density
    """
    tau_0 = (2*J+1) * const.h * const.c**3 * A * N / (16 * np.pi * J * const.k_B * nu**2 * dv * Tex) \
            * (np.exp(const.h*nu/(const.k_B*Tex) - 1)
               * np.exp(-(J+1)*const.h*nu/(2*const.k_B*Tex)))

    return tau_0.to(1)

def tau_0_ratio_f(J1, J2, A1, A2, nu1, nu2, Tex):
    """
    Calculates the optical depth ratio of two transitions 1 and 2
    """
    tau0_ratio = J2 * (2*J1+1) * A1 / (J1 * (2*J2+1) * A2) \
                 * (nu2 / nu1)**2 \
                 * (np.exp(const.h*nu1/(const.k_B*Tex)) - 1) / (np.exp(const.h*nu2/(const.k_B*Tex)) - 1) \
                 * np.exp(const.h * ((J2+1)*nu2 - (J1+1)*nu1) / (2*const.k_B*Tex))
    return tau0_ratio


freq_caract_CO = {'1-0':115.27120180*u.GHz,
                  '2-1':230.53800000*u.GHz,
                  '3-2':345.79598990*u.GHz}
def Inu_f(TR, nu):
    """
    Intensity as a function of radiation temperature TR
    """
    Inu = 2 * const.k_B * nu**2 * TR / const.c**2
    return Inu.to(u.Jy) / u.sr

def Jybeam2K(nu, Inu, maja, mina):
    """
    Calculates the intensity in K given the beamsize major axis and minor axis
    of the beam and the intensity in fluxunit/beam
    Parameters
    ----------
    nu: astropy.units.Unit
        Frequency of the transition
    Inu: astropy.units.Unit
        Intensity, can be given in u.Jy/beam or u.mJy/beam
    maja: astropy.units.Unit
        Major axis of the beam (for example in u.arcsec)
    mina: astropy.units.Unit
        Minor axis of the beam (for example in u.arcsec)
    """
    beamsize_sr = mb_sa_gaussian_f(maja, mina)
    return Inu_2_T0(nu, Inu/beamsize_sr)

def K2Jybeam(nu, TR, maja, mina):
    """
    Calculates the intensity in Jy/beam given the beamsize major axis and minor axis
    of the beam and the intensity in Kelvin
    Parameters
    ----------
    nu: astropy.units.Unit
        Frequency of the transition
    TR: astropy.units.Unit
        Intensity in Kelvin
    maja: astropy.units.Unit
        Major axis of the beam (for example in u.arcsec)
    mina: astropy.units.Unit
        Minor axis of the beam (for example in u.arcsec)
    """
    beamsize_sr = mb_sa_gaussian_f(maja, mina)
    return Inu_f(TR, nu) * beamsize_sr

def intensityK2Tb(nu, Jnu):
    """
    Calculates the brightness temperature given the intensity in temperature
    nu: astropy.units.Unit
        Frequency of the transition
    Jnu: astropy.units.Unit
       Intensity in temperature units (K)
    """
    Tb = const.h * nu / (const.k_B * np.log(const.h*nu/(const.k_B*Jnu)+1))
    return Tb.to(u.K)

def pixel2offset(pixels, cdelts, center_pixels, bounds=False):
    """
    Calculates the offset coordinates of pixels with respect to center_pixels.
    If bounds=True it returns the extent in offsets, thought to use it in
    imshow.
    """
    if bounds:
        if (pixels[0]==center_pixels[0]) or (pixels[1]==center_pixels[1]):
            print("Pixel to convert to offset should not be the center pixel when bounds=True!")
            return None
        else:
            xaddpix = 0 if pixels[0]<center_pixels[0] else 1
            yaddpix = 0 if pixels[1]<center_pixels[1] else 1
            offset = [
                cdelts[0]*(pixels[0]+xaddpix-center_pixels[0]-0.5),
                cdelts[1]*(pixels[1]+yaddpix-center_pixels[1]-0.5),
            ]
    else:
         offset = [
            cdelts[0]*(pixels[0]-center_pixels[0]),
            cdelts[1]*(pixels[1]-center_pixels[1]),
        ]
    return offset
