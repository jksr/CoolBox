from typing import Union

from collections import OrderedDict

import svgutils.compose as sc
import matplotlib.pyplot as plt

from coolbox.utilities.figtools import cm2inch
from coolbox.utilities.filetool import get_uniq_tmp_file
from coolbox.utilities import GenomeRange
from coolbox.core.track import Track
from coolbox.core.frame import Frame
from .base import SuperFrame



from coolbox.core.track.arcs.fetch import FetchParix
from coolbox.utilities.bed import process_bedpe
from matplotlib.path import Path as mpath
from matplotlib.patches import PathPatch
import pandas as pd




class CoLinearView(SuperFrame):
    """Compose two track and a matching track.

    Parameters
    ----------
    middle : Track
        Center track for show 'contact map-like' plot,
        should support for plotting images with 2d GenomeRanges.

    top : Frame, optional
        Frame plot in the top of the center track.


    bottom : Frame, optional

    middle_width : float
        The width of center track, unit in cm. default 20.

    space : float
        Space between frame and center, unit in cm. default 0.5
    """

    def __init__(self,
                 middle: Track,
                 top: Union[Frame, Track] = None,
                 bottom: Union[Frame, Track] = None,
                 **kwargs,
                 ):
        sub_frames = OrderedDict({
            "top": top,
            "bottom": bottom,
        })
        sub_frames = {k: v for k, v in sub_frames.items() if v is not None}
        self.__check_sub_frames(middle, sub_frames)
        # explicitely use matrix style
#         center.properties['style'] = 'matrix'
        properties = {
            "sub_frames": sub_frames,
            "middle_track": middle,
            "middle_height": 5,
            "middle_width": 20,
#             "trbl": "1212",
            "space": 1,
            "cm2px": 28.5,
            "padding_left": 1,
        }
        properties.update(**kwargs)

        super().__init__(properties)
        self.__adjust_sub_frames_width()

    def cm2px(self, vec):
        return [i * self.properties['cm2px'] for i in vec]

    @staticmethod
    def __check_sub_frames(middle, sub_frames):
#         from ..frame import Frame
#         from ...track.base import Track
        sub_f_names = ", ".join(sub_frames.keys())
#       TODO
#         if (not isinstance(middle, Track)) and (not hasattr(center, "plot_joint")):
#             raise TypeError("center should be a Track type instance with plot_joint method, "
#                             "for example Cool, DotHiC, ...")
        for k, f in sub_frames.items():
            if not isinstance(f, Frame):
                if isinstance(f, Track):
                    sub_frames[k] = (Frame() + f)  # convert track to frame
                else:
                    raise TypeError(f"{sub_f_names} should be Frame object.")

    def __adjust_sub_frames_width(self):
        for k, f in self.properties['sub_frames'].items():
            width_ratios = f.properties['width_ratios']
            middle_ratio = width_ratios[1]
            new_width = self.properties['middle_width'] / middle_ratio
            f.properties['width'] = new_width
    def plot_middle(self, gr1: GenomeRange, gr2: GenomeRange):
        middle_track = self.properties['middle_track']
        width = cm2inch(self.properties['middle_width'])
        height = cm2inch(self.properties['middle_height'])
        fig, ax = plt.subplots(figsize=(width, height))
        
        middle_track.plot(ax, gr1, gr2)
        middle_track.plot_coverages(ax, gr1, gr2)
        ax.set_axis_off()
        path = get_uniq_tmp_file(prefix='center', suffix='.svg')
#         fig.subplots_adjust(wspace=0, hspace=0.0, left=0, right=1, bottom=0, top=1)
        fig.subplots_adjust(wspace=0, hspace=0.0, left=0, right=1, bottom=0, top=1)
        fig.savefig(path)
        plt.close()
        return sc.SVG(path)

    def frame_granges(self, gr1=None, gr2=None):
        self.goto(gr1, gr2)
        gr1, gr2 = self.current_range
        return {
            #'top':(gr1,gr1),
            #'bottom':(gr2,gr2),
            'top':(gr1,None),
            'bottom':(gr2,None),
        }

    def plot(self, gr1=None, gr2=None):
        """

        Parameters
        ----------
        gr1 : {str, GenomeRange}
            First genome range

        gr2 : {str, GenomeRange}, optional
            Second genome range
        """
        frame2grange = self.frame_granges(gr1, gr2)
        gr1, gr2 = self.current_range
        sub_frames = self.properties['sub_frames']

        frame_svgs = self.plot_frames(frame2grange)
        center_svg = self.plot_middle(gr1, gr2)

        center_offsets = self.__get_center_offsets(sub_frames)

        center_svg.move(*self.cm2px(center_offsets))
        self.__transform_sub_svgs(frame_svgs, sub_frames, center_offsets)

        figsize = self.cm2px(self.__get_figsize(sub_frames))
        return sc.Figure(f"{figsize[0]}px", f"{figsize[1]}px",
                        sc.Panel(center_svg),
                        *[sc.Panel(svg) for svg in frame_svgs.values()])

    def fetch_data(self, gr1=None, gr2=None) -> dict:
        """

        Parameters
        ----------
        gr1 : {str, GenomeRange}, optional
            First genome range

        gr2 : {str, GenomeRange}, optional
            Second genome range
        """
        frame2grange = self.frame_granges(gr1, gr2)
        gr1, gr2 = self.current_range
        sub_frames = self.properties['sub_frames']

        tracks_data = {
            pos: sub_frames[pos].fetch_data(
                frame2grange[pos][0], gr2=frame2grange[pos][1]
            )
            for pos, fr in sub_frames.items()
        }

        tracks_data['middle'] = self.properties['middle_track'].fetch_data(gr1, gr2=gr2)

        return tracks_data

    def __transform_sub_svgs(self, sub_svgs, sub_frames, center_offsets):
        c_width = self.properties['middle_width']
        c_height = self.properties['middle_height']
        space = self.properties['space']
        pd_left = self.properties['padding_left']
        if 'top' in sub_svgs:
            s = sub_svgs['top']
            offsets = [pd_left, 0]
            if 'left' in sub_svgs:
                offsets[0] += sub_frames['left'].properties['height'] + space
            s.move(*self.cm2px(offsets))
        if 'bottom' in sub_svgs:
            s = sub_svgs['bottom']
            offsets = [pd_left, c_height]
            if 'left' in sub_svgs:
                offsets[0] += sub_frames['left'].properties['height'] + space
            if 'top' in sub_svgs:
                offsets[1] += sub_frames['top'].properties['height']
            offsets = self.cm2px(offsets)
            s.move(*offsets)


    def __get_center_offsets(self, sub_frames):
        space = self.properties['space']
        center_offsets = [self.properties['padding_left'], 0]  # x, y (left, top)
        if 'top' in sub_frames:
            f = sub_frames['top']
            center_offsets[0] += f.properties['width'] * f.properties['width_ratios'][0]
            center_offsets[1] += space + f.properties['height']
        if 'bottom' in sub_frames:
            f = sub_frames['bottom']
            if 'top' not in sub_frames:
                center_offsets[0] += f.properties['width'] * f.properties['width_ratios'][0]
        return center_offsets

    def __get_figsize(self, sub_frames):
        space = self.properties['space']
        middle_width = self.properties['middle_width']
        middle_height = self.properties['middle_height']
        size = [middle_width, middle_height]  # width, height
        if 'top' in sub_frames:
            f = sub_frames['top']
            size[0] = f.properties['width']
            size[1] += f.properties['height'] + space
        if 'bottom' in sub_frames:
            f = sub_frames['bottom']
            size[0] = f.properties['width']
            size[1] += f.properties['height'] + space
        self.properties['width'] = size[0] #+ self.properties['padding_left']
        self.properties['height'] = size[1]
        return size

    def add_track(self, track, pos=None):
        pass
#         print('hhhh')
#         sub_frames = self.properties['sub_frames']
#         if pos is None:
#             pos = list(sub_frames.keys())[-1]
#         frame = sub_frames[pos]
#         frame.add_track(track)

class CoXAxis(Track, FetchParix):
    """
    The x axis track.
    Parameters
    ----------
    height : float, optional
        Height of Spacer track. (Default: XAxis.DEFAULT_HEIGHT)
    fontsize : int, optional
        Font size of XAxis. (Default: XAxis.DEFAULT_FONTSIZE)
    name (str, optional):
        Track's name.
    """

    #TODO style options
    DEFAULT_PROPERTIES = {
        "height": 5,
        "fontsize": 15,
        "patchpath_style":dict(linestyle=":", 
                               linewidth=1,
                               edgecolor="#bdbdbd",
                               facecolor="#9ecae1",
                               alpha=0.5,),
    }

    FIELDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2",
              "name", "score", "strand1", "strand2"]    

    def __init__(self, file, **kwargs):
        properties = CoXAxis.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(properties)
        self.bgz_file = process_bedpe(file)

    def fetch_data(self, gr: GenomeRange, gr2: GenomeRange, **kwargs) -> pd.DataFrame:
        # filter peaks manually for hicpeaks style in fetch_plot_data
        df = self.fetch_intervals(self.bgz_file, gr, gr2)
        # TODO the returned df has no named columns, may cause error
        if len(df) == 0:
            return df

        columns = list(df.columns)
        for i, col in enumerate(self.FIELDS):
            if i >= len(columns):
                break
            columns[i] = col
        df.columns = columns
        for col in ['start1', 'end1', 'start2', 'end2']:
            df[col] = df[col].astype(int)
        return df
    
    def fetch_plot_data(self, gr: GenomeRange, gr2: GenomeRange, **kwargs) -> pd.DataFrame:
        df = self.fetch_data(gr,gr2, **kwargs)
        if len(df)!=0:
            df['name'] = df['strand1']==df['strand2']
        return df[['chrom1','start1','end1',
                   'chrom2','start2','end2','name']].rename(columns={'name':'same direction'})
    
    @staticmethod
    def _genomerange_to_ticklabels(ticks,gr):
        a = (gr.end-gr.start)/(ticks[-1]-ticks[0])
        b = gr.start-a*ticks[0]
        labels = [a*x+b for x in ticks]
        
        if labels[-1] - labels[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in labels]
            labels[-2] += " Kb"

        elif 1e5 < labels[-1] - labels[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in labels]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in labels]
            labels[-2] += " Mbp"
            
        return labels

    @staticmethod
    def _ticks_to_ticklabels(labels):
        if labels[-1] - labels[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in labels]
            labels[-2] += " Kb"

        elif 1e5 < labels[-1] - labels[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in labels]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in labels]
            labels[-2] += " Mbp"
            
        return labels
    
    @staticmethod
    def _genomerange_to_axisrange(gr,refgr,ticks):
        a = (ticks[-1]-ticks[0])/(refgr.end-refgr.start)
        b = ticks[0]-a*refgr.start
        return (a*x+b for x in gr)
    
    
    def plot(self, ax, gr1: GenomeRange, gr2: GenomeRange, **kwargs):
        self.ax = ax
        
        ax1 = ax.secondary_xaxis('top',)
        ax1.tick_params(axis='x', which='major', labelsize=int(self.properties['fontsize']))
        ax1.tick_params(axis='x', which='minor', top='on')
        
        ax2 = ax.secondary_xaxis('bottom',)
        ax2.tick_params(axis='x', which='major', labelsize=int(self.properties['fontsize']))
        ax2.tick_params(axis='x', which='minor', bottom='on')
        
        
        ticks = ax.get_xticks()
        ticks1 = ax1.get_xticks()
        ticks2 = ax2.get_xticks()

        labels1 = CoXAxis._genomerange_to_ticklabels(ticks1,gr1)
        labels2 = CoXAxis._genomerange_to_ticklabels(ticks2,gr2)
        
        ax1.set_xticklabels(labels1)
        ax2.set_xticklabels(labels2)
        
        
#         ax1.set_axis_off()
#         ax2.set_axis_off()
#         ax.margins(y=0)
        
        df = self.fetch_plot_data(gr1, gr2,**kwargs)
        

        base1=1
        base2=0

        for _,row in df.iterrows():
            _,s1,e1,_,s2,e2,same_d = row
            s1,e1 = CoXAxis._genomerange_to_axisrange((s1,e1),gr1,ticks)
            s2,e2 = CoXAxis._genomerange_to_axisrange((s2,e2),gr2,ticks)

            if same_d:
                pathcode, pathverts = zip(*[(mpath.MOVETO,(s1,base1)),
                                            (mpath.LINETO,(e1,base1)),
                                            (mpath.LINETO,(e2,base2)),
                                            (mpath.LINETO,(s2,base2)),
                                            (mpath.CLOSEPOLY,(s1,base1))])
            else:
                pathcode, pathverts = zip(*[(mpath.MOVETO,(s1,base1)),
                                            (mpath.LINETO,(e1,base1)),
                                            (mpath.LINETO,(s2,base2)),
                                            (mpath.LINETO,(e2,base2)),
                                            (mpath.CLOSEPOLY,(s1,base1))])
            ax.add_patch(PathPatch(mpath(pathverts, pathcode ),
                                   **self.properties["patchpath_style"]))
            
    def plot_bad(self, ax, gr1: GenomeRange, gr2: GenomeRange, **kwargs):
        import numpy as np
        self.ax = ax
        xlim = ax.get_xlim()

        self.ax.set_xlim(gr1.start, gr1.end)
        labels1 = ax.get_xticks()
        self.ax.set_xlim(gr2.start, gr2.end)
        labels2 = ax.get_xticks()

        self.ax.set_xlim(*xlim)
        
        
        ax1 = ax.secondary_xaxis('top',)
        ax1.tick_params(axis='x', which='major', labelsize=int(self.properties['fontsize']))
        ax1.tick_params(axis='x', which='minor', top='on')
        
        ax2 = ax.secondary_xaxis('bottom',)
        ax2.tick_params(axis='x', which='major', labelsize=int(self.properties['fontsize']))
        ax2.tick_params(axis='x', which='minor', bottom='on')
        
        
        ticks = ax.get_xticks()
        
        ticks1 = np.linspace(0,1,len(labels1))
        ax1.set_xticks(ticks1)
        ticks2 = np.linspace(0,1,len(labels2))
        ax2.set_xticks(ticks2)

        labels1 = CoXAxis._ticks_to_ticklabels(labels1)
        labels2 = CoXAxis._ticks_to_ticklabels(labels2)
        
        ax1.set_xticklabels(labels1)
        ax2.set_xticklabels(labels2)
        
        
#         ax1.set_axis_off()
#         ax2.set_axis_off()
#         ax.margins(y=0)
        
        df = self.fetch_plot_data(gr1, gr2,**kwargs)
        

        base1=1
        base2=0

        for _,row in df.iterrows():
            _,s1,e1,_,s2,e2,same_d = row
            s1,e1 = CoXAxis._genomerange_to_axisrange((s1,e1),gr1,ticks)
            s2,e2 = CoXAxis._genomerange_to_axisrange((s2,e2),gr2,ticks)

            if same_d:
                pathcode, pathverts = zip(*[(mpath.MOVETO,(s1,base1)),
                                            (mpath.LINETO,(e1,base1)),
                                            (mpath.LINETO,(e2,base2)),
                                            (mpath.LINETO,(s2,base2)),
                                            (mpath.CLOSEPOLY,(s1,base1))])
            else:
                pathcode, pathverts = zip(*[(mpath.MOVETO,(s1,base1)),
                                            (mpath.LINETO,(e1,base1)),
                                            (mpath.LINETO,(s2,base2)),
                                            (mpath.LINETO,(e2,base2)),
                                            (mpath.CLOSEPOLY,(s1,base1))])
            ax.add_patch(PathPatch(mpath(pathverts, pathcode ),
                                   **self.properties["patchpath_style"]))
