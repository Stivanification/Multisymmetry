using Gadfly
using Compose
using DataFrames

using Colors
base03=parse(Colorant,"#002b36");
base02=parse(Colorant,"#073642");
base01=parse(Colorant,"#586e75");
base00=parse(Colorant,"#657b83");
base0=parse(Colorant,"#839496");
base1=parse(Colorant,"#839496");
base2=parse(Colorant,"#eee8d5");
base3=parse(Colorant,"#fdf6e3");

yellow=parse(Colorant,"#b58900");
orange=parse(Colorant,"#cb4b16");
red=parse(Colorant,"#dc322f");
magenta=parse(Colorant,"#d33682");
violet=parse(Colorant,"#6c71c4");
blue=parse(Colorant,"#268bd2");
cyan=parse(Colorant,"#3aa198");
green=parse(Colorant,"#859900");

dash = 6 * Compose.cm
dot = 1.5 * Compose.cm
gap = 1 * Compose.cm

x = rand(10)

l1=layer(x=x,y=[1:10;],Geom.line,Theme(line_style=[dot]))
l2=layer(x=x,y=ones(x),Geom.line,Theme(line_style=[dash]))
l3=layer(x=x,y=2*ones(x),Geom.line,Theme(line_style=[gap]))
l4=layer(x=x,y=3*ones(x),Geom.line,Theme(line_style=[dot,dash,dot]))
l5=layer(x=x,y=4*ones(x),Geom.line,Theme(line_style=[gap,dash]))


plot(dataset(x=[1:10;],y=[1:10;],
    Theme(highlight_width=0pt, # lose the white ring around the points
    default_point_size=3pt,    # size of the dots
    default_color=green,     # color of the dots
    background_color=base02,   # ... background color ...
    grid_color=base3,     # the lines
    minor_label_color=base3,  # numbers on axes color
    key_label_color=base3, # color key numbers color
    key_title_color=cyan,  # label of color key color
    major_label_color=cyan, # axes label and title color
    major_label_font_size=20pt, # font size for axes label and title
    minor_label_font_size=15pt, #font size for numbers on axes
	panel_opacity = 1
);)