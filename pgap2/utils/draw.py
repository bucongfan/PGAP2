import os
import tempfile
import numpy as np
from loguru import logger
from itertools import chain
from collections import Counter, OrderedDict

import pyecharts.options as opts
from pyecharts.globals import ThemeType
from pyecharts.commons.utils import JsCode

from pgap2.lib.species import Species
from pgap2.utils.supply import run_command
from decimal import Decimal, getcontext


preprocess_cfg_dict = [{"cid": "PGAP2_preprocess_ani_scatter3d", "width": "518px", "height": "258px", "top": "67px", "left": "586px"},
                       {"cid": "PGAP2_preprocess_species_pie", "width": "531px",
                           "height": "259px", "top": "65px", "left": "39px"},
                       {"cid": "PGAP2_preprocess_gene_completeness_bar", "width": "532px",
                           "height": "257px", "top": "343px", "left": "38px"},
                       {"cid": "PGAP2_preprocess_gene_length_box", "width": "530px",
                           "height": "259px", "top": "613px", "left": "38px"},
                       {"cid": "PGAP2_preprocess_genome_content_bar", "width": "520px",
                           "height": "256px", "top": "343px", "left": "587px"},
                       {"cid": "PGAP2_preprocess_gene_code_heatmap", "width": "500px",
                           "height": "530px", "top": "344px", "left": "1133px"},
                       {"cid": "PGAP2_preprocess_half_core_line", "width": "523px", "height": "258px", "top": "613px", "left": "585px"}]

postprocess_cfg_dict = [{"cid": "PGAP2_postprocess_group_freq_line", "width": "685px", "height": "408px", "top": "33.5px", "left": "53px"},
                        {"cid": "PGAP2_postprocess_stat_pan_scatter", "width": "683px",
                            "height": "408px", "top": "35px", "left": "759px"},
                        {"cid": "PGAP2_postprocess_mean_line", "width": "687px",
                            "height": "366px", "top": "864.5px", "left": "47px"},
                        {"cid": "PGAP2_postprocess_uni_line", "width": "687px",
                            "height": "370px", "top": "1258px", "left": "46px"},
                        {"cid": "PGAP2_postprocess_min_line", "width": "691px",
                            "height": "366px", "top": "865.5px", "left": "757px"},
                        {"cid": "PGAP2_postprocess_var_line", "width": "693px",
                            "height": "371px", "top": "1258px", "left": "757px"},
                        {"cid": "PGAP2_postprocess_profile_box", "width": "687px",
                            "height": "362px", "top": "464.5px", "left": "51px"},
                        {"cid": "PGAP2_postprocess_new_clusters_box", "width": "685px", "height": "365px", "top": "464px", "left": "757px"}]


def get_cord(Z, strain_name, strain_order):
    for i, name in enumerate(strain_order):
        if name == strain_name:
            return Z[i]
    raise ValueError(f"Cannot find {strain_name} in ani dict")


def preprocess_species_pie(darb_dict, sp):
    from pyecharts.charts import Pie
    inner_data = []
    outer_data = []
    num_outlier_ani = len(darb_dict['ani'])
    num_outlier_sgl = len(darb_dict['single_gene'])
    num_strain = len(sp.strain_dict)
    num_pass = num_strain-num_outlier_ani-num_outlier_sgl
    inner_data.append(['Total', num_strain])
    outer_data.append(['pass', num_pass])
    outer_data.append(['outlier-ani', num_outlier_ani])
    outer_data.append(['outlier-single', num_outlier_sgl])
    pie = Pie()
    pie.add(
        series_name="Species",
        data_pair=inner_data,
        radius=[0, "30%"],
        label_opts=opts.LabelOpts(position="inner"),
    ).add(series_name="Filtration",
          radius=["40%", "55%"],
          data_pair=outer_data,
          label_opts=opts.LabelOpts(
              position="outside",
              formatter="{a|{a}}{abg|}\n{hr|}\n {b|{b}: }{c}  {per|{d}%}  ",
              background_color="#eee",
              border_color="#aaa",
              border_width=0,
              border_radius=4,
              rich={
                  "a": {"color": "#999", "lineHeight": 22, "align": "center"},
                  "abg": {
                      "backgroundColor": "#e3e3e3",
                      "width": "100%",
                      "align": "right",
                      "height": 22,
                      "borderRadius": [4, 4, 0, 0],
                  },
                  "hr": {
                      "borderColor": "#aaa",
                      "width": "100%",
                      "borderWidth": 0.5,
                      "height": 0,
                  },
                  "b": {"fontSize": 16, "lineHeight": 33},
                  "per": {
                      "color": "#eee",
                      "backgroundColor": "#334455",
                      "padding": [2, 4],
                      "borderRadius": 2,
                  },
              },
          ),
          ).set_series_opts(
        tooltip_opts=opts.TooltipOpts(
            trigger="item", formatter="{a} <br/>{b}: {c} ({d}%)"
        )
    )
    pie.set_global_opts(title_opts=opts.TitleOpts(title="Filtration result"),
                        legend_opts=opts.LegendOpts(is_show=False),
                        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name':
                                                                               'pgap2.preprocessing.SpeciesStatisticsPlot'}}))
    return pie
    pie.render('pie.html')
    exit()


def preprocess_ani_graph(outlier_dict, sp: Species):
    from pyecharts.charts import Graph

    ani_dict = sp.ani_dict
    nodes = []
    links = []
    categories = [{'name': 'Pass'},
                  {'name': 'Gene count outlier'},
                  {'name': 'ANI outlier'},
                  {'name': 'Both outlier'}]
    darb_strain = sp.get_darb()
    darb_strain_name = sp.strain_dict[darb_strain].strain_name
    nodes.append(opts.GraphNode(name=darb_strain_name,
                 symbol_size=10, category='Pass', value=1))
    for strain in ani_dict:
        if strain == darb_strain:
            continue
        if strain in outlier_dict['single_gene']:
            category = 'Gene count outlier'
        elif strain in outlier_dict['ani']:
            category = 'ANI outlier'
        else:
            category = 'Pass'
        if strain in outlier_dict['single_gene'] and strain in outlier_dict['ani']:
            category = 'Both outlier'
        strain_name = sp.strain_dict[strain].strain_name
        nodes.append(opts.GraphNode(name=strain_name,
                     symbol_size=10, category=category, value=ani_dict[strain]))
        links.append(opts.GraphLink(source=strain_name,
                     target=darb_strain_name))

    graph = Graph(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    graph.add("", nodes, links, categories, repulsion=50,
              edge_length=50, linestyle_opts=opts.LineStyleOpts(curve=0.2))
    graph.set_global_opts(title_opts=opts.TitleOpts(title="ANI relationship"),
                          legend_opts=opts.LegendOpts(is_show=False),
                          toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right", feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']}, "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name': 'pgap2.preprocessing.AniGraph'}}))
    return graph
    graph.render('graph.html')

# def preprocess_ani_scatter3d(Z, darb_dict, ani_dict):
#     from pyecharts.charts import Scatter3D
#     scatter3D = Scatter3D(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
#     strain_order = list(sorted(ani_dict.keys()))
#     darb_score_list = []
#     for darb in darb_dict:
#         this_z = []  # cord of this species
#         species = darb_dict[darb]['species']
#         for strain in darb_dict[darb]['strain']:
#             darb_score = darb_dict[darb]['darb_score'][strain]
#             this_cord = get_cord(Z, strain, strain_order)
#             this_cord = this_cord.tolist()
#             is_outlier = 1 if strain in darb_dict[darb]['outlier'] else 0
#             is_darb = 1 if strain == darb else 0
#             print(f"{strain} {is_outlier} {this_cord}")
#             # the darb_score used to map the size of the circle
#             darb_score_list.append(darb_score)
#             this_cord.append(round(darb_score, 4))
#             this_cord.append(strain)  # strain name
#             this_cord.append(is_outlier)  # is the outlier of the species
#             this_cord.append(is_darb)  # is the darb of the species
#             this_z.append(this_cord)
#         ###############
#         # darb_dict[darb]['strain']里的strain不包含darb本身，因此需要单独增加darb的坐标
#         darb_score = darb_dict[darb]['darb_score'][darb]
#         this_cord = get_cord(Z, darb, strain_order)
#         this_cord = this_cord.tolist()
#         is_outlier = 1 if strain in darb_dict[darb]['outlier'] else 0
#         is_darb = 1
#         darb_score_list.append(darb_score)
#         this_cord.append(round(darb_score, 4))
#         this_cord.append(darb)  # strain name
#         this_cord.append(is_outlier)  # is the outlier of the species
#         this_cord.append(is_darb)  # is the darb of the species
#         this_z.append(this_cord)
#         ################
#         scatter3D.add(series_name=species, data=this_z, grid3d_opts=opts.Grid3DOpts(
#             width=100, height=100, depth=100, splitarea_opts=opts.SplitAreaOpts(is_show=True)),)
#     scatter3D.set_global_opts(visualmap_opts=[opts.VisualMapOpts(is_show=False, type_="size", dimension=3, min_=min(darb_score_list), max_=max(darb_score_list))],
#                               title_opts=opts.TitleOpts(
#                                   title="ANI relationship",),
#                               legend_opts=opts.LegendOpts(
#                                   pos_top="10%", pos_left="2%"), toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right", feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']}, "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name': 'pgap2.preprocessing.StrainAniDistancePlot'}}),
#                               tooltip_opts=opts.TooltipOpts(is_show=True, trigger_on="mousemove|click", axis_pointer_type='cross', formatter=JsCode(
#                                   """function(params){return 'x: '+params.data[0]+'<br/>'+'y: '+params.data[1]+'<br/>'+'z: '+params.data[2]+'<br/>'+'darb score: '+params.data[3]+'<br/>'+'strain: '+params.data[4]+'<br/>'+'is_outlier: '+params.data[5]+'<br/>'+'is_darb: '+params.data[6]}"""
#                               )))

#     return scatter3D
#     scatter3D.render("scatter3d.html")


def calculate_edges(min_val, max_val, num_bins):
    edges = np.linspace(min_val, max_val, num_bins + 1).tolist()
    return sorted([round(_, 5) for _ in edges])


def calculate_logedges(min_val, max_val, num_bins):
    edges = np.logspace(np.log10(min_val), np.log10(
        max_val), num_bins + 1).tolist()
    return sorted(edges)
    # edges = np.linspace(min_val, max_val, num_bins + 1).tolist()
    # return sorted([round(_, 5) for _ in edges])


def calculate_counts(data, edges):
    # 初始化计数字典
    result = OrderedDict({edge: 1 for edge in edges[:-1]})
    for value in data:
        for i in range(len(edges) - 1):
            if edges[i] <= value < edges[i + 1]:
                result[edges[i]] += 1
                break
             # 最后一个区间改为左闭右闭
            if value == edges[-1]:
                result[edges[-2]] += 1
    return list(result.values())


def edges_for_xaxis(edges):
    xaxis = []
    for i in range(len(edges)-1):
        if i == len(edges)-2:
            xaxis.append(f'[{edges[i]},{edges[i+1]}]')
        else:
            xaxis.append(f'[{edges[i]},{edges[i+1]})')
    return xaxis


def postprocess_attr_line(benchmark_stat: dict, para_id, attr_name, ofh):
    from pyecharts.charts import Line
    line = Line(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    real_edges = []
    if attr_name in ('min', 'mean'):
        edges = calculate_edges(para_id, 1, 10)
        real_edges = edges_for_xaxis(edges)
    elif attr_name == 'var':
        pseudo_min = 1
        for pangroup, pangroup_stat in benchmark_stat.items():
            for i, value in enumerate(pangroup_stat):
                if value < pseudo_min and value != 0:
                    pseudo_min = value
        pseudo_min = pseudo_min/10
        # substitute 0 with a pseudo value in benchmark_stat
        for pangroup, pangroup_stat in benchmark_stat.items():
            for i, value in enumerate(pangroup_stat):
                if value == 0:
                    pangroup_stat[i] = pseudo_min
        edges = calculate_logedges(
            pseudo_min, max(chain.from_iterable(benchmark_stat.values())), 10)
        real_edges = edges_for_xaxis(edges)

    elif attr_name == 'uni':
        edges = calculate_edges(para_id, 1, 10)
        real_edges = [f'<{para_id}'] + edges_for_xaxis(edges)
        edges.insert(0, 0)

    line.add_xaxis(real_edges)
    for pangroup, pangroup_stat in benchmark_stat.items():
        if pangroup_stat:
            counts = calculate_counts(pangroup_stat, edges)
            for edge, count in zip(edges, counts):
                if attr_name == 'var':
                    edge = f"{edge:.3e}"
                ofh.write(f"{attr_name}\t{pangroup}\t{edge}\t{count}\n")
            line.add_yaxis(pangroup, counts,
                           linestyle_opts=opts.LineStyleOpts(width=2),)
            #   category_gap=0, stack='stack1')
            logger.debug(f"{attr_name} {pangroup} {counts}")

    line.set_series_opts(
        label_opts=opts.LabelOpts(
            is_show=False,
        )
    )

    line.set_global_opts(tooltip_opts=opts.TooltipOpts(
        is_show=True, trigger="axis", axis_pointer_type="cross",
    ),
        xaxis_opts=opts.AxisOpts(
            type_="category",
            name_rotate=45,
            axispointer_opts=opts.AxisPointerOpts(
                is_show=True, type_="shadow"),
    ),
        yaxis_opts=opts.AxisOpts(
            is_show=True, axispointer_opts=opts.AxisPointerOpts(is_show=False), is_scale=True, type_='log'),
        title_opts=opts.TitleOpts(title=attr_name,),
        legend_opts=opts.LegendOpts(pos_left="40%", border_width=0),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right=0, feature={"dataView": {'title': 'view the data',
                                                                                                                   'lang': ['Data view', 'Close', 'Refresh']},
                                                                                                      "saveAsImage": {'type_': 'png', 'title': 'save as png',
                                                                                                                      'pixel_ratio': 5, 'name': 'pgap2.postprocess.uniPlot'}}),
    )
    return line
    line.render('test_attr.html')


def postprocess_mci_bar(benchmark_stat, iter_list):
    from pyecharts.charts import Bar
    lable = sorted(iter_list, reverse=True)

    bar = Bar(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    bar.add_xaxis(lable)
    for type in benchmark_stat:
        counts = Counter(benchmark_stat[type])
        stat_list = []
        for _ in lable:
            if _ not in counts:
                stat_list.append(0)
            else:
                stat_list.append(counts[_])

        bar.add_yaxis(type, stat_list,
                      stack='stack1', category_gap="50%")

    bar.set_series_opts(
        label_opts=opts.LabelOpts(
            is_show=False,
        )
    )

    bar.set_global_opts(tooltip_opts=opts.TooltipOpts(
        is_show=True, trigger="axis", axis_pointer_type="cross",
    ),
        xaxis_opts=opts.AxisOpts(
            type_="category",
            axispointer_opts=opts.AxisPointerOpts(
                is_show=True, type_="shadow"),
    ),
        yaxis_opts=opts.AxisOpts(
            is_show=True, axispointer_opts=opts.AxisPointerOpts(is_show=False),),
        title_opts=opts.TitleOpts(title="MCI",),
        legend_opts=opts.LegendOpts(pos_left="40%", border_width=0), toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right=0, feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']}, "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name': 'pgap2.postprocess.mciPlot'}}),
    )
    return bar
    bar.render('mci_bar.html')


def preprocess_benchmark_bar(benchmark_stat, outdir):
    from pyecharts.charts import Bar
    lable = sorted(benchmark_stat['total'].keys(), reverse=True)
    total = [benchmark_stat['total'][_] for _ in lable]
    benchmark = [benchmark_stat['benchmark'][_] for _ in lable]
    paralog = [benchmark_stat['paralog'][_] for _ in lable]
    bar = Bar(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    bar.add_xaxis(lable)
    bar.add_yaxis("total_clust", total)
    bar.add_yaxis("bench_clust", benchmark)
    bar.add_yaxis("para_clust", paralog)

    bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False))

    bar.set_global_opts(tooltip_opts=opts.TooltipOpts(
        is_show=True, trigger="axis", axis_pointer_type="cross"
    ),
        xaxis_opts=opts.AxisOpts(
            type_="category",
            axispointer_opts=opts.AxisPointerOpts(
                is_show=True, type_="shadow"),
    ), datazoom_opts=[opts.DataZoomOpts(type_="inside", xaxis_index=[0, 1])],
        title_opts=opts.TitleOpts(title="Cluster component",),
        legend_opts=opts.LegendOpts(pos_left="40%"), toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right", feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']}, "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name': 'pgap2.preprocessing.ClusterComponentPlot'}}),
    )
    bar.render(f'{outdir}/cluster_component.html')
    return bar


def preprocess_half_core_line(benchmark_stat, sp):
    from pyecharts.charts import Line
    strain_name_list = []
    half_core = []
    single_cloud = []
    for strain_index in range(len(benchmark_stat)):
        strain_name = sp.strain_dict[strain_index].strain_name
        strain_name_list.append(strain_name)
        half_core.append(benchmark_stat[strain_index]['half_core'])
        single_cloud.append(benchmark_stat[strain_index]['single_cloud'])

    line = Line()
    line.add_xaxis(xaxis_data=strain_name_list)
    line.add_yaxis(
        series_name="Half core",
        y_axis=half_core,
        markpoint_opts=opts.MarkPointOpts(
            data=[
                opts.MarkPointItem(type_="max", name="Maximal value"),
                opts.MarkPointItem(type_="min", name="Minimal value"),
            ]
        ),
        markline_opts=opts.MarkLineOpts(
            data=[opts.MarkLineItem(type_="average", name="Average value")]
        ),
    )
    line.add_yaxis(
        series_name="Single cloud",
        y_axis=single_cloud,
        markpoint_opts=opts.MarkPointOpts(
            data=[
                opts.MarkPointItem(type_="max", name="Maximal value"),
                opts.MarkPointItem(type_="min", name="Minimal value"),
            ]
        ),
        markline_opts=opts.MarkLineOpts(
            data=[opts.MarkLineItem(type_="average", name="Average value")]
        ),
    )
    line.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    line.set_global_opts(
        title_opts=opts.TitleOpts(title="Gene cluster"),
        tooltip_opts=opts.TooltipOpts(trigger="axis"),
        xaxis_opts=opts.AxisOpts(type_="category", boundary_gap=False),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                               'name': 'pgap2.preprocessing.HalfOrSingleGeneClusterPlot'}}),
    )
    return line
    line.render("half_count.html")


def preprocess_gene_completeness_bar(genome_attr_dict, sp):
    from pyecharts.charts import Bar
    complete_list = []
    incomplete_list = []
    species_list = []
    for strain in genome_attr_dict:
        complete = genome_attr_dict[strain]['total_gene_num']
        incomplete = genome_attr_dict[strain]['value_incomplete']
        total = complete+incomplete
        complete_list.append(
            {"value": complete, "percent": complete/total})
        incomplete_list.append(
            {"value": incomplete, "percent": incomplete/total})
        strain_name = sp.strain_dict[strain].strain_name
        species_list.append(strain_name)
    bar = Bar(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    bar.add_xaxis(species_list)
    bar.add_yaxis("complete",
                  complete_list, stack="stack1", category_gap="50%")
    bar.add_yaxis("incomplete", incomplete_list,
                  stack="stack1", category_gap="50%")
    bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    # bar.set_series_opts(
    #     label_opts=opts.LabelOpts(
    #         position="right",
    #         formatter=JsCode(
    #             "function(x){return Number(x.data.percent * 100).toFixed() + '%';}"
    #         ),
    #     ))
    bar.set_global_opts(tooltip_opts=opts.TooltipOpts(
        is_show=True, trigger="axis", axis_pointer_type="cross"
    ),
        xaxis_opts=opts.AxisOpts(
            type_="category",
            axispointer_opts=opts.AxisPointerOpts(
                is_show=True, type_="shadow"),
    ),
        yaxis_opts=opts.AxisOpts(
            name="#gene",
            type_="value",
            axislabel_opts=opts.LabelOpts(formatter="{value}"),
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),
    ), datazoom_opts=[opts.DataZoomOpts(type_="inside", xaxis_index=[0, 1])],
        title_opts=opts.TitleOpts(title="Gene completeness",),
        legend_opts=opts.LegendOpts(pos_left="40%"),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                               'name': 'pgap2.preprocessing.GenomeCompletenessPlot'}}),
    )

    return bar
    bar.render('gene_completeness.html')
    exit()


def preprocess_gene_length_box(genome_attr_dict, sp: Species):
    from pyecharts.charts import Boxplot, Line
    score_list = []
    species_list = []
    co_list = []
    for strain_index in sp.strain_dict:
        this_list = sp.gene_len[strain_index]
        co_list.append(np.mean(this_list))
        # sample_count = len(this_list) if len(this_list) < 300 else 300
        # this_list = random.sample(this_list, sample_count)
        score_list.append(this_list)
        strain_name = sp.strain_dict[strain_index].strain_name
        species_list.append(strain_name)

    # co_list = [round(_/max(co_list), 2) for _ in co_list]
    box = Boxplot()
    box.add_xaxis(xaxis_data=species_list)
    box.add_yaxis(
        series_name="gene length",
        y_axis=box.prepare_data(score_list),
        tooltip_opts=opts.TooltipOpts(
            formatter=JsCode(
                """function(param) { return [
                            'Strain ' + param.name + ': ',
                            'upper: ' + param.data[0],
                            'Q1: ' + param.data[1],
                            'median: ' + param.data[2],
                            'Q3: ' + param.data[3],
                            'lower: ' + param.data[4]
                        ].join('<br/>') }"""
            )
        ),
    )
    box.set_global_opts(
        title_opts=opts.TitleOpts(title="Average gene length"),
        legend_opts=opts.LegendOpts(pos_left="40%"),
        tooltip_opts=opts.TooltipOpts(
            trigger="item", axis_pointer_type="shadow"),
        xaxis_opts=opts.AxisOpts(
            name_gap=30,
            boundary_gap=True,
            splitarea_opts=opts.SplitAreaOpts(
                areastyle_opts=opts.AreaStyleOpts(opacity=1)
            ),
            axislabel_opts=opts.LabelOpts(formatter="{value}"),
            splitline_opts=opts.SplitLineOpts(is_show=False),
        ),
        yaxis_opts=opts.AxisOpts(
            name="#base",
            type_="value",
            axislabel_opts=opts.LabelOpts(formatter="{value}"),
            splitline_opts=opts.SplitLineOpts(is_show=True),
        ),
        datazoom_opts=[
            opts.DataZoomOpts(type_="inside", xaxis_index=0),
        ], toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                         feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                                  "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                                  'name': 'pgap2.preprocessing.GeneLength'}})
    )

    line = Line().add_xaxis(xaxis_data=species_list).add_yaxis(
        series_name="Average gene length",
        y_axis=co_list,
        label_opts=opts.LabelOpts(is_show=False),
    )
    box.overlap(line)

    return box
    box.render('gene_length.html')


def preprocess_genome_content_bar(genome_attr_dict, sp):
    from pyecharts.charts import Bar
    species_list = []
    n_list = []
    a_list = []
    t_list = []
    c_list = []
    g_list = []

    for strain_index in sp.strain_dict:
        content_dict = genome_attr_dict[strain_index]['content']
        a_list.append(content_dict['a'])
        t_list.append(content_dict['t'])
        c_list.append(content_dict['c'])
        g_list.append(content_dict['g'])
        n_list.append(content_dict['n'])
        strain_name = sp.strain_dict[strain_index].strain_name
        species_list.append(strain_name)

    bar = Bar(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    bar.add_xaxis(species_list)
    bar.add_yaxis("A",
                  a_list, stack="stack1", category_gap="50%", )
    bar.add_yaxis("T",
                  t_list, stack="stack1", category_gap="50%", )
    bar.add_yaxis("C",
                  c_list, stack="stack1", category_gap="50%", )
    bar.add_yaxis("G",
                  g_list, stack="stack1", category_gap="50%", )
    bar.add_yaxis("N",
                  n_list, stack="stack1", category_gap="50%", )

    bar.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    bar.set_global_opts(toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                                               'name': 'pgap2.preprocessing.GenomeCompositionPlot'}}),
                        tooltip_opts=opts.TooltipOpts(
        is_show=True, trigger="axis", axis_pointer_type="cross"
    ),
        xaxis_opts=opts.AxisOpts(
            type_="category",
            axispointer_opts=opts.AxisPointerOpts(
                is_show=True, type_="shadow"),
    ),
        yaxis_opts=opts.AxisOpts(
            name="#base",
            type_="value",
            axislabel_opts=opts.LabelOpts(formatter="{value}"),
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),
    ), datazoom_opts=[opts.DataZoomOpts(type_="inside")],
        title_opts=opts.TitleOpts(title="Genome content"), legend_opts=opts.LegendOpts(pos_left="40%"))

    return bar
    bar.render('genome_content.html')


def preprocess_gene_code_heatmap(sp):
    from pyecharts.charts import HeatMap

    species_list = []
    heatmap_list = []
    max_coden = 0

    for y, code in enumerate(sp.gene_code):
        for strain_index, count in sp.gene_code[code].items():
            heatmap_list.append([strain_index, y, count])
            max_coden = max(max_coden, count)
    ref_code = list(sp.gene_code.keys())
    species_list = [
        sp.strain_dict[_].strain_name for _ in range(len(sp.strain_dict))]

    heat = HeatMap().add_xaxis(species_list)
    heat.add_yaxis(
        "",
        ref_code,
        heatmap_list,
        label_opts=opts.LabelOpts(is_show=False, position="inside"),
    )
    # heat.set_series_opts(label_opts=opts.LabelOpts(is_show=True))
    heat.set_global_opts(title_opts=opts.TitleOpts(title="Gene coden distribution"),
                         legend_opts=opts.LegendOpts(pos_left="15%"), visualmap_opts=opts.VisualMapOpts(is_show=True, orient="horizontal", pos_top="1%", pos_right="5%", max_=max_coden),
                         datazoom_opts=[opts.DataZoomOpts(
                             type_="inside"), ], toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                                                               feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                                                                        "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                                                                        'name': 'pgap2.preprocessing.StrainCodenUsagePlot'}})
                         )
    return heat
    heat.render("gene_code_heatmap.html")


def postprocess_newclusters_box(pan_profile):
    '''
    pan_profile=[[sampling number],[core_list],[pan_list]]
    '''
    from pyecharts.charts import Boxplot
    boxplot = Boxplot()
    boxplot.add_xaxis(xaxis_data=pan_profile[0])
    boxplot.add_yaxis("", y_axis=boxplot.prepare_data(pan_profile[1]))
    boxplot.set_global_opts(
        title_opts=opts.TitleOpts(title="New clusters"),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                               'name': 'pgap2.postprocessing.NewClustersPlot'}}),
        xaxis_opts=opts.AxisOpts(
            name_location='middle', name_gap=30,
            name_textstyle_opts=opts.TextStyleOpts(font_size='14'),
            name='Strain number', type_="category",
            axispointer_opts=opts.AxisPointerOpts(is_show=True, type_="shadow"),),
        yaxis_opts=opts.AxisOpts(
            name='New Gene Cluster Number',
            type_="value",
            name_location='middle', name_gap=50,
            min_=min(pan_profile[1][-1]),
            name_textstyle_opts=opts.TextStyleOpts(font_size='14',),
            axistick_opts=opts.AxisTickOpts(is_show=False),
            splitline_opts=opts.SplitLineOpts(is_show=False),),
        # tooltip_opts=opts.TooltipOpts(is_show=True, trigger="axis", axis_pointer_type="cross"),
        legend_opts=opts.LegendOpts(),
    )

    return boxplot
    boxplot.render('pangenome_newclusters.boxplot.html')


def postprocess_profile_box(pan_profile):
    '''
    pan_profile=[[sampling number],[core_list],[pan_list]]
    '''
    from pyecharts.charts import Boxplot
    boxplot = Boxplot().add_xaxis(xaxis_data=pan_profile[0])
    boxplot.add_yaxis(
        "Core genome", y_axis=boxplot.prepare_data(pan_profile[1]))
    boxplot.add_yaxis(
        "Pan genome", y_axis=boxplot.prepare_data(pan_profile[2]))
    boxplot.set_global_opts(
        title_opts=opts.TitleOpts(title="Pangenome curve"),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']},
                                               "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5,
                                                               'name': 'pgap2.postprocessing.PanGenomeCurvePlot'}}),
        xaxis_opts=opts.AxisOpts(
            name_location='middle', name_gap=30,
            name_textstyle_opts=opts.TextStyleOpts(font_size='14'),
            name='Strain number', type_="category",
            axispointer_opts=opts.AxisPointerOpts(is_show=True, type_="shadow"),),
        yaxis_opts=opts.AxisOpts(
            name='Clusters',
            type_="value",
            name_location='middle', name_gap=50,
            min_=min(pan_profile[1][-1]),
            name_textstyle_opts=opts.TextStyleOpts(font_size='14',),
            axistick_opts=opts.AxisTickOpts(is_show=False),
            splitline_opts=opts.SplitLineOpts(is_show=False),),
        # tooltip_opts=opts.TooltipOpts(is_show=True, trigger="axis", axis_pointer_type="cross"),
        legend_opts=opts.LegendOpts(),
    )

    return boxplot
    boxplot.render('pangenome_profile.boxplot.html')


def postprocess_group_freq(group_freq):
    from pyecharts.charts import Line, Pie

    '''
    group_freq=[[group],[freq], {pangroup}]
    '''
    pan_group = group_freq[2]
    lineplot = Line()
    lineplot.add_xaxis([str(f) for f in group_freq[0]],)
    # lineplot.add_yaxis(
    #     '', group_freq[1], category_gap=0, color='#73c0df')
    lineplot.add_yaxis(
        '', group_freq[1], is_smooth=True, color='#73c0df', areastyle_opts=opts.AreaStyleOpts(opacity=1, color='#73c0de'),)

    lineplot.set_series_opts(

        label_opts=opts.LabelOpts(is_show=False),
        markarea_opts=opts.MarkAreaOpts(
            is_silent=False,
            data=[
                opts.MarkAreaItem(name="Cloud", x=(
                    "0.0", "0.15"), label_opts=opts.LabelOpts(is_show=False), itemstyle_opts=opts.ItemStyleOpts(color="red", opacity=0.1),),
                opts.MarkAreaItem(name="Shell", x=(
                    "0.15", "0.95"), label_opts=opts.LabelOpts(is_show=False), itemstyle_opts=opts.ItemStyleOpts(color="yellow", opacity=0.1),),
                opts.MarkAreaItem(name="Core", x=(
                    "0.95", "1.0"), label_opts=opts.LabelOpts(is_show=False), itemstyle_opts=opts.ItemStyleOpts(color="blue", opacity=0.1),),
            ])
    )
    lineplot.set_global_opts(
        title_opts=opts.TitleOpts(title="Strain frequency"),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right",
                                      feature={"dataView": {'title': 'view the data', 'lang': ['Data view', 'Close', 'Refresh']}, "saveAsImage": {'type_': 'png', 'title': 'save as png', 'pixel_ratio': 5, 'name': 'pgap2.postprocessing.StrainFrequencyPlot'}}),
        xaxis_opts=opts.AxisOpts(name_location='middle', name_gap=30, name_textstyle_opts=opts.TextStyleOpts(font_size='14'),
                                 name='Strain frequency',
                                 boundary_gap=False,
                                 type_="category",
                                 axispointer_opts=opts.AxisPointerOpts(is_show=True, type_="shadow"),),
        yaxis_opts=opts.AxisOpts(
            name='Clusters',
            type_="value",
            name_location='middle', name_gap=50,
            name_textstyle_opts=opts.TextStyleOpts(font_size='14',),
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),),
        axispointer_opts=opts.AxisPointerOpts(is_show=False),
        tooltip_opts=opts.TooltipOpts(
            is_show=True, trigger="axis", axis_pointer_type="cross"),
        legend_opts=opts.LegendOpts(
            border_width=0, pos_top="4%", is_show=False),
    )
    pie = (
        Pie()
        # .set_colors(['#5470c6', '#91cc75', '#fac858', '#ee6666'])
        .add(
            series_name="pan-group",
            data_pair=[
                ['Strict core [1.0, 1.0]', pan_group['Strict_core']],
                ['Core [0.99, 1.0)', pan_group['Core']],
                ['Soft core [0.95, 0.99)', pan_group['Soft_core']],
                ['Shell (0.15, 0.95)', pan_group['Shell']],
                ['Cloud [0.0, 0.15]', pan_group['Cloud']],
            ],
            center=["50%", "45%"],
            radius="38%",
        )
        .set_series_opts(tooltip_opts=opts.TooltipOpts(is_show=True, trigger="item", formatter="{a} <br/>{b}: {c} ({d}%)"))
    )

    lineplot.overlap(pie)
    return lineplot
    lineplot.render('pangenome_groupfreq.lineplot.html')


def postprocess_stat_para(para_dict: dict):
    '''
    Draw postprocess plot

    :param para_dict: dict, {pangroup: [[x],[y]]}

    return:
    html obj
    '''
    from pyecharts.charts import Scatter, Line
    scatter = Scatter()
    scatter.add_dataset
    for pangroup, coord in para_dict.items():  # prefix in (gene, strain, expect)
        x = coord[0]
        y = coord[1]
        scatter.add_xaxis(xaxis_data=x,)
        scatter.add_yaxis(pangroup, y_axis=y,
                          label_opts=opts.LabelOpts(is_show=False),
                          )

    scatter.set_global_opts(
        title_opts=opts.TitleOpts(title="Paralogs distribution"),
        legend_opts=opts.LegendOpts(pos_left="40%", border_width=0),
        toolbox_opts=opts.ToolboxOpts(orient='horizontal', pos_bottom="bottom", pos_right="right", feature={"dataView": {'title': 'view the data', 'lang': [
                                      'Data view', 'Close', 'Refresh']},
            "saveAsImage": {'type_': 'png', 'title': 'save as png',
                            'pixel_ratio': 5, 'name': 'pgap2.postprocessing.ParaDistributionPlot'}}),
        tooltip_opts=opts.TooltipOpts(
            is_show=True, trigger="axis", axis_pointer_type="cross"),
        xaxis_opts=opts.AxisOpts(
            name='Paralogous strain', name_location='middle', name_gap=30,
            name_textstyle_opts=opts.TextStyleOpts(font_size='14',), type_="value",
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),
        ),

        yaxis_opts=opts.AxisOpts(
            name='Paralogous gene', name_location='middle', name_gap=30,
            name_textstyle_opts=opts.TextStyleOpts(font_size='14',), type_="value",
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),),
    )

    # line=Line()
    # line.add_xaxis(xaxis_data=[1,max(x)])
    # line.add_yaxis('One paralog per strain',[1,max(x)],label_opts=opts.LabelOpts(is_show=False))

    return scatter
    scatter.render('postprocess_pan_dinuq.html')


# def postprocess_pan_dinuq(data_dict: dict):
#     from pyecharts.charts import Scatter
#     scatter = Scatter()
#     scatter.add_dataset
#     # this_df=self.df.loc[:,['per_pan','conformity','best_k','this_discore','this_cmass_x','this_cmass_y']]
#     max_strain_num = 0
#     for prefix in data_dict:
#         x_list = []
#         y_list = []
#         discore_list = []
#         for index, pan in enumerate(data_dict[prefix]):
#             pan_name = pan['pan_clust_name']
#             strain_num = pan['strain_num']
#             gene_num = pan['gene_num']
#             max_strain_num = strain_num if strain_num > max_strain_num else max_strain_num
#             x_list.append(float(pan['cmass'][0]))
#             this_y = []
#             this_y.append(float(pan['cmass'][1]))
#             this_y.append(pan_name)
#             this_y.append(strain_num)
#             this_y.append(gene_num)
#             this_y.append(float(pan['radius']))
#             this_y.append(str(pan['clust_from']))
#             this_y.append(float(pan['best_k']))
#             this_y.append(float(pan['prob']))
#             this_y.append(float(pan['discore']))
#             this_y.append(float(pan['conformity']))

#             discore_list.append(float(pan['discore']))

#             y_list.append(this_y)
#         scatter.add_xaxis(xaxis_data=x_list,)
#         scatter.add_yaxis(prefix, y_axis=y_list,
#                           label_opts=opts.LabelOpts(is_show=False))
#     scatter.set_global_opts(toolbox_opts=opts.ToolboxOpts(orient='horizontal',pos_bottom="bottom",pos_right="right",feature={"dataView": {'title':'view the data','lang':['Data view', 'Close', 'Refresh']},"saveAsImage":{'type_':'png','title':'save as png','pixel_ratio':5,'name':'pgap2.postprocessing.DiscoreDistributionPlot'}}),tooltip_opts=opts.TooltipOpts(
#         formatter=JsCode(
#             "function (params) {return 'Pan: '+params.value[2] +'<br/>'+'Strain_num: '+params.value[3] +'<br/>'+'Gene_num: '+params.value[4] +'<br/>'+'Probability: '+params.value[8] +'<br/>'+ 'Radius: '+params.value[5] +'<br/>' +'From_cluster: '+params.value[6] +'<br/>'+'Discore: '+params.value[9] +'<br/>'+ 'Conformity: '+params.value[10];}"
#         )
#     ),
#         visualmap_opts=[opts.VisualMapOpts(is_show=True, type_="color", dimension=4, min_=0, max_=max_strain_num, orient="horizontal", pos_top='1%', pos_right='10%'), opts.VisualMapOpts(
#             is_show=False, type_="size", dimension=5, max_=max(discore_list), min_=min(discore_list))],
#         xaxis_opts=opts.AxisOpts(
#             name='Vector 1', name_location='middle', name_gap=30, name_textstyle_opts=opts.TextStyleOpts(font_size='14',), type_="value", splitline_opts=opts.SplitLineOpts(is_show=True)
#     ),
#         yaxis_opts=opts.AxisOpts(
#             name='Vector 2',
#             type_="value",
#             name_location='middle', name_gap=30, name_textstyle_opts=opts.TextStyleOpts(font_size='14',),
#             axistick_opts=opts.AxisTickOpts(is_show=True),
#             splitline_opts=opts.SplitLineOpts(is_show=True),
#     ), legend_opts=opts.LegendOpts(border_width=0))
#     return scatter
#     scatter.render('postprocess_pan_dinuq.html')

def preprocess_draw_vector(**kwargs):
    outdir = kwargs['outdir']
    single_file = kwargs['single_file']
    ani_threshold = kwargs['ani_threshold']
    sfw = kwargs['sfw']
    input_prep = os.path.join(outdir, 'preprocess.stat.tsv')
    input_code = os.path.join(outdir, 'preprocess.gene_code.csv')
    if not os.path.exists(input_prep):
        logger.error(f'No preprocess.stat.tsv found in {outdir}')
        raise FileNotFoundError(
            f'No preprocess.stat.tsv found in {outdir}')
    if not os.path.exists(input_code):
        logger.error(f'No preprocess.gene_code.csv found in {outdir}')
        raise FileNotFoundError(
            f'No preprocess.gene_code.csv found in {outdir}')
    single_file = '--single_file' if single_file else ''
    run_command(
        f"Rscript {sfw} -a {input_prep} -b {input_code} -o {outdir} --ani_thre {ani_threshold} {single_file}")


def postprocess_draw_vector(**kwargs):
    target = kwargs['target']
    outdir = kwargs['outdir']
    single_file = kwargs['single_file']
    sfw = kwargs['sfw']
    if target == 'stat':
        input_file = os.path.join(outdir, 'postprocess.stat_attrs.tsv')
        if not os.path.exists(input_file):
            logger.error(f'No postprocess.stat_attrs.tsv found in {outdir}')
            raise FileNotFoundError(
                f'No postprocess.stat_attrs.tsv found in {outdir}')
        single_file = '--single_file' if single_file else ''
        run_command(f"Rscript {sfw} -a {input_file} -o {outdir} {single_file}")
    elif target == 'profile':
        input_files = (os.path.join(outdir, 'postprocess.pan_group_stat.tsv'),
                       os.path.join(
                           outdir, 'postprocess.clust_strain_freq.tsv'),
                       os.path.join(
            outdir, 'postprocess.rarefaction.tsv'),
            os.path.join(
            outdir, 'postprocess.new_clusters.tsv'),
            os.path.join(outdir, 'postprocess.para_stat.tsv'))
        for input_file in input_files:
            if not os.path.exists(input_file):
                logger.error(f'No {input_file} found in {outdir}')
                raise FileNotFoundError(f'No {input_file} found in {outdir}')
        (input_a, input_b, input_c, input_d, input_e) = input_files
        single_file = '--single_file' if single_file else ''
        run_command(
            f"Rscript {sfw} -a {input_a} -b {input_b} -c {input_c} -d {input_d} -e {input_e} -o {outdir} {single_file}")
    else:
        logger.error(f'Invalid target: {target}')


def postprocess_draw(**kwargs):
    '''
    Draw postprocess plot

    return:
    html file
    '''
    target = kwargs['target']

    from pyecharts.charts import Page
    page = Page(layout=Page.DraggablePageLayout)

    if target == 'stat':
        basic = kwargs.get('basic')
        outdir = kwargs.get('outdir')
        group_freq = kwargs.get('group_freq')
        pan_para_stat = kwargs.get('pan_para_stat')
        pan_profile = kwargs.get('pan_profile')
        new_clusters = kwargs.get('new_clusters')

        pan_profile_box = postprocess_profile_box(pan_profile=pan_profile)
        new_clusters_box = postprocess_newclusters_box(new_clusters)
        group_freq_line = postprocess_group_freq(group_freq)
        stat_pan_scatter = postprocess_stat_para(pan_para_stat)

        ofh = open(f"{outdir}/postprocess.stat_attrs.tsv", 'w')
        ofh.write(f"Attr\tGroup\tEdge\tCount\n")
        var_line = postprocess_attr_line(
            basic.var_dict, basic.para_id, 'var', ofh)
        mean_line = postprocess_attr_line(
            basic.mean_dict, basic.para_id, 'mean', ofh)
        min_line = postprocess_attr_line(
            basic.min_dict, basic.para_id, 'min', ofh)
        uni_line = postprocess_attr_line(
            basic.uni_dict, basic.para_id, 'uni', ofh)
        ofh.close()

        pan_profile_box.chart_id = 'PGAP2_postprocess_profile_box'
        new_clusters_box.chart_id = 'PGAP2_postprocess_new_clusters_box'
        group_freq_line.chart_id = 'PGAP2_postprocess_group_freq_line'
        stat_pan_scatter.chart_id = 'PGAP2_postprocess_stat_pan_scatter'
        mean_line.chart_id = 'PGAP2_postprocess_mean_line'
        uni_line.chart_id = 'PGAP2_postprocess_uni_line'
        min_line.chart_id = 'PGAP2_postprocess_min_line'
        var_line.chart_id = 'PGAP2_postprocess_var_line'

        page.add(group_freq_line)
        page.add(stat_pan_scatter)
        page.add(mean_line)
        page.add(uni_line)
        page.add(min_line)
        page.add(var_line)
        page.add(pan_profile_box)
        page.add(new_clusters_box)
        page.page_title = 'PGAP2 Postprocess Stat'

    elif target == 'profile':
        group_freq = kwargs.get('group_freq')
        pan_para_stat = kwargs.get('pan_para_stat')
        pan_profile = kwargs.get('pan_profile')
        new_clusters = kwargs.get('new_clusters')
        outdir = kwargs.get('outdir')

        pan_profile_box = postprocess_profile_box(pan_profile=pan_profile)
        new_clusters_box = postprocess_newclusters_box(new_clusters)
        group_freq_line = postprocess_group_freq(group_freq)
        stat_pan_scatter = postprocess_stat_para(pan_para_stat)

        group_freq_line.chart_id = 'PGAP2_postprocess_group_freq_line'
        stat_pan_scatter.chart_id = 'PGAP2_postprocess_stat_pan_scatter'
        pan_profile_box.chart_id = 'PGAP2_postprocess_profile_box'
        new_clusters_box.chart_id = 'PGAP2_postprocess_new_clusters_box'

        page.add(group_freq_line)
        page.add(stat_pan_scatter)
        page.add(pan_profile_box)
        page.add(new_clusters_box)
        page.page_title = 'PGAP2 Postprocess Profile'

    output_html = f"{outdir}/postprocess_{target}.html"

    fd, fname = tempfile.mkstemp(dir=outdir, suffix='.html')
    try:
        page.render(fname)
        Page.save_resize_html(
            fname, cfg_dict=postprocess_cfg_dict, dest=output_html)
    finally:
        os.remove(fname)
    return output_html


def preprocess_draw(outlier_dict, outdir, sp, genome_attr_dict):
    from pyecharts.charts import Page
    # ani_scatter3d = preprocess_ani_scatter3d(darb_dict, sp.ani_dict)
    species_pie = preprocess_species_pie(outlier_dict, sp)
    ani_graph = preprocess_ani_graph(outlier_dict, sp)
    gene_completeness_bar = preprocess_gene_completeness_bar(
        genome_attr_dict, sp)
    half_core_line = preprocess_half_core_line(genome_attr_dict, sp)
    gene_length_box = preprocess_gene_length_box(genome_attr_dict, sp)
    genome_content_bar = preprocess_genome_content_bar(genome_attr_dict, sp)
    gene_code_heatmap = preprocess_gene_code_heatmap(sp)

    ani_graph.chart_id = 'PGAP2_preprocess_ani_scatter3d'
    species_pie.chart_id = 'PGAP2_preprocess_species_pie'
    gene_completeness_bar.chart_id = 'PGAP2_preprocess_gene_completeness_bar'
    gene_length_box.chart_id = 'PGAP2_preprocess_gene_length_box'
    genome_content_bar.chart_id = 'PGAP2_preprocess_genome_content_bar'
    gene_code_heatmap.chart_id = 'PGAP2_preprocess_gene_code_heatmap'
    half_core_line.chart_id = 'PGAP2_preprocess_half_core_line'

    page = Page(layout=Page.DraggablePageLayout)
    page.add(ani_graph)
    page.add(species_pie)
    page.add(gene_completeness_bar)
    page.add(gene_length_box)
    page.add(genome_content_bar)
    page.add(gene_code_heatmap)
    page.add(half_core_line)
    page.page_title = 'PGAP2 Preprocess'

    fd, fname = tempfile.mkstemp(dir=outdir, suffix='.html')
    try:
        page.render(fname)
        Page.save_resize_html(
            fname, cfg_dict=preprocess_cfg_dict, dest=f"{outdir}/preprocess_results.html")
    finally:
        os.remove(fname)
    return f"{outdir}/pgap2.preprocess.html"
