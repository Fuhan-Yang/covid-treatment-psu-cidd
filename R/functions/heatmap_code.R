### Preamble
suppressWarnings({
  suppressMessages({
    require(ggplot2)
    require(tidyr)
    require(gridExtra)
    require(RColorBrewer)
    require(plyr)
    require(dplyr)
    require(cowplot)
    require(reshape)
    require(reshape2)
    require(grid)
    require(scales)
    require(metR)
  })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_model_sims <- args[1]
input_analysis_tools <- args[2]
input_viz_tools <- args[3]
input_multipanel_flag <-args[4] # if "FALSE", do not make multipanel fig
input_arrow_flag <- args[5] # if "FALSE", no arrows on plot, if "TRUE" add arrows to plot
input_other_levels_flag <- args[6] # if "FALSE" do not create plots for other levels, if "TRUE" create and save
output_path <- args[7]



source(input_viz_tools)
source(input_analysis_tools)
out = read.csv(input_model_sims)
out[which(round(out$Ds,8) == round(1/24,8)),"Ds"] = 1/24 # adjust for truncation when saving
#contours = read.csv(input_contours)

input_multipanel_flag = as.logical(input_multipanel_flag)
input_arrow_flag = as.logical(input_arrow_flag)
input_other_levels_flag = as.logical(input_other_levels_flag)

#### contours ####
# calculate contours
out_sub = out %>% filter(objective == "tot.cases",tot.tests>0, Ts>0)
contours = data.frame(tot.tests = integer(),
                      sens = double(),
                      level = integer(),
                      x = double(),
                      y = double())
for(i in unique(out_sub$d)){
  for(j in unique(out_sub$sens)){
    c = contour(out_sub %>% filter(d == i, sens == j), "Ds", "tot.tests", "median", c(250,500,1000))
    if(!is.null(c)){contours = rbind(contours,cbind(i,j,c))}
  }
}
colnames(contours) = c("d", "sens","level", "Ds", "tot.tests")


#### heatmap - cases ####
if(input_multipanel_flag){
  ## Panel A,B
  out_sub =out %>% filter(objective == "tot.cases",sens == 1,tot.tests>0, Ts>0)
  y_range = c(min((out_sub[, "median"])),max(out_sub[ "median"]))
  
  out_plot = out %>% filter(objective == "tot.cases",tot.tests>0, Ts>0)
  
  
  p1=ggplot(out_plot %>% filter(tot.tests %in% c(100,500,2000,5000), Ds != 1/24, sens == 1))+
    geom_raster(aes(y = d, x= Ds, fill = median))+
    geom_contour(aes(y = d, x = Ds, z = median), breaks = c(250,500,1000,2000), color ="black")+
    geom_text_contour(aes(y = d, x = Ds, z = median), stroke = 0.1, breaks = c(250,500,1000,2000))+
    facet_grid(cols = vars(tot.tests), labeller = labeller(tot.tests = labs_tot_tests))+
    scale_fill_distiller(palette = "YlOrBr",name = "Cumulative infections", direction = 0, limits = y_range, labels = comma)+
    scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)")+
    scale_y_continuous(expand = c(0,0), name = "NPI intensity", labels = percent_format(accuracy = 1))+
    theme_classic()+
    theme(legend.position = "bottom",
          legend.key.width = unit(1.2,"cm"),
          panel.spacing = unit(1.2, "lines"))
  
  p2=ggplot(out_plot %>% filter(Ds %in% c(1/24, 1,2,5), tot.tests %in% seq(500,5000,500), sens == 1))+
    geom_tile(aes(y = d, x= tot.tests/10000, fill = median))+
    geom_contour(aes(y = d, x = tot.tests/10000, z = median), breaks = c(250,500,1000,2000), color ="black")+
    geom_text_contour(aes(y = d, x = tot.tests/10000, z = median), stroke = 0.1, breaks = c(250,500,1000,2000))+
    facet_grid(cols = vars(Ds), labeller = labeller(Ds = labs_Ds))+
    scale_fill_distiller(palette = "YlOrBr",name = "Cumulative infections", direction = 0, limits = y_range)+
    scale_x_continuous(expand = c(0,0), name = "Test administration\n (% of population tested daily)", labels = percent_format(accuracy = 1))+
    scale_y_continuous(expand = c(0,0), name = "NPI intensity", labels = percent_format(accuracy = 1))+
    theme_classic()+
    theme(legend.position = "bottom",
          legend.key.width = unit(1.2,"cm"),
          panel.spacing = unit(1.2, "lines"))
  
  ## Panel C
  labels = data.frame(labels = paste0(seq(0,30,10),"%"),
                      y = 0.51,
                      contours %>% filter(tot.tests == 5000, level == 500, round(d,4) %in% round(seq(0,0.3,0.1),4)) %>% 
                        group_by(d) %>% 
                        summarise(Ds = mean(Ds)))
  
  labels2 = data.frame(x = min(labels$Ds)-1,
                       y =  min(labels$y),
                       lab = "NPI: ")
  
  p = ggplot(out %>% filter(objective == "tot.cases", Ts>0, round(d,2) %in% round(seq(0,0.3,0.1),2)),
             aes(x = Ds, y = tot.tests/10000))+
    geom_contour(aes(z = median, color = d, linetype = as.factor(sens), group = paste(d, sens)), breaks = c(500), size = 1.2)+
    geom_text(data = labels, aes(x = Ds, y = y, label = labels, color = d))+
    geom_text(data = labels2, aes(x = x, y =y, label = lab), color = "black", hjust = 0)+
    guides(color = FALSE, lty = guide_legend(override.aes = list(col = 'black')))+
    scale_color_distiller(palette = "Reds", direction = 0, name = "NPI Intensity", labels = percent_format(accuracy = 1), limits = c(0,0.55))+
    scale_linetype_manual(values = c("dashed", "solid"), labels = c("90%","100%"), name = "Effective\nSensitivity")+
    scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)", limits = c(0,8))+
    scale_y_continuous(expand = c(0,0), name = "Test administration\n (% of population tested daily)", labels = percent_format(accuracy = 1), limits = c(0,0.52))+
    theme_classic()+
    theme(legend.direction = "vertical",
          legend.key.width = unit(1.2,"cm"),
          panel.spacing = unit(1.2, "lines"))
  
  
  if(input_arrow_flag == "TRUE"){
    
    contours$round_d = round(contours$d, 4)
    contours_plot = contours %>%
      filter(round_d %in% c(0,0.1,0.2,0.3), level ==500) 
    
    start_sens = 1
    start_d = 0.3
    start_tot_tests = 1400
    start_df = contours %>% filter(sens == start_sens, round_d == start_d, level == 500)
    start = as.data.frame(t(unlist(approx(start_df$tot.tests, start_df$Ds, start_tot_tests))))
    colnames(start) = c("tot.tests","Ds")
    start$sens= start_sens
    start$d = start_d
    
    # lower sens, same d
    sub = contours_plot %>% filter(sens == start_sens-0.1, round_d == start_d)
    contour1 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = start$Ds))))
    contour2 = as.data.frame(t(unlist(approx(sub$tot.tests,sub$Ds,xout = start$tot.tests))))
    contour_tmp = rbind(contour1, c(contour2$y, contour2$x))
    contour3 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = min(contour_tmp$x[1],contour_tmp$x[2]) + abs(contour_tmp$x[1] - contour_tmp$x[2])/2))))
    contour_tmp = rbind(contour_tmp, contour3)
    contour_tmp$thickness = "a"
    contour_thick = cbind(start, contour_tmp)
    
    # same sens, lower d
    sub = contours_plot %>% filter(sens == start_sens, round_d == round(start_d-0.1,4))
    contour1 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = start$Ds))))
    contour2 = as.data.frame(t(unlist(approx(sub$tot.tests,sub$Ds,xout = start$tot.tests))))
    contour_tmp = rbind(contour1, c(contour2$y, contour2$x))
    contour3 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = min(contour_tmp$x[1],contour_tmp$x[2]) + abs(contour_tmp$x[1] - contour_tmp$x[2])/2))))
    contour_tmp = rbind(contour_tmp, contour3)
    contour_tmp$thickness = "b"
    #contour = rbind(contour,contour_tmp)
    contour_medium = cbind(start, contour_tmp)
    
    # lower sens, lower d
    sub = contours_plot %>% filter(sens == start_sens-0.1, round_d == round(start_d-0.1,1))
    contour1 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = start$Ds))))
    contour2 = as.data.frame(t(unlist(approx(sub$tot.tests,sub$Ds,xout = start$tot.tests))))
    contour_tmp = rbind(contour1, c(contour2$y, contour2$x))
    contour3 = as.data.frame(t(unlist(approx(sub$Ds, sub$tot.tests,xout = min(contour_tmp$x[1],contour_tmp$x[2]) + abs(contour_tmp$x[1] - contour_tmp$x[2])/2))))
    contour_tmp = rbind(contour_tmp, contour3)
    contour_tmp$thickness = "c"
    #contour = rbind(contour,contour_tmp)
    contour_light = cbind(start, contour_tmp)
    
    # for text
    # contour_medium[1,"y"] - contour_medium[1,"tot.tests"]
    # contour_medium[2,"Ds"] - contour_medium[2,"x"]
    # 
    # contour_thick[1,"y"] - contour_thick[1,"tot.tests"]
    # contour_thick[2,"Ds"] - contour_thick[2,"x"]
    # (contour_thick[2,"Ds"] - contour_thick[2,"x"]) * 24
    
    p = p + 
      geom_point(data = as.data.frame(start), aes(x = Ds, y= tot.tests/10000), color = "grey55", size = 3.5)+
      geom_segment(data = contour_thick, aes(xend = x, yend = y/10000), 
                   color = "grey55", lineend = "round", linejoin = "round",
                   size = 1.5, arrow = arrow(length = unit(0.25,"cm")))+
      geom_segment(data = contour_medium, aes(xend = x, yend = y/10000), 
                   color = "grey55", lineend = "round", linejoin = "round",
                   size = 1.05, arrow = arrow(length = unit(0.225,"cm")))+
      geom_segment(data = contour_light, aes(xend = x, yend = y/10000), 
                   color = "grey55", lineend = "round", linejoin = "round",
                   size = 0.7, arrow = arrow(length = unit(0.225,"cm")))
    
  }
  
  l = get_legend(p)
  p = p + 
    theme(legend.position = "none")+
    annotation_custom(grob=l,ymin = 0, ymax=0.2, xmin=5, xmax=8)
  
  
  #### make full multipanel figure ####
  
  # combine panels a and b, add legend
  l_col = get_legend(p1)
  f2 = plot_grid(p2+theme(legend.position = "none"), 
                 p1+theme(legend.position = "none"), 
                 align = "vh",labels = c("A", "B"),ncol =1)
  f2 = plot_grid(f2, l_col, ncol  =1, rel_heights = c(0.9,0.1))
  
  # add panel c
  f = plot_grid(f2, plot_grid(p, ggplot()+theme_void(),ncol =1, rel_heights = c(0.9,0.1)), 
                ncol = 2, labels = c(NA,"C"),rel_widths = c(0.6,0.4))
  ggsave(file.path(output_path),f,width = 15, height = 7, units = "in")
}






#### plot for other levels ####

if(input_other_levels_flag){
  contours$round_d = round(contours$d, 4)
  contours_plot = contours %>%
    filter(round_d %in% c(0,0.1,0.2,0.3), level %in% c(250,1000)) 
  
  text = data.frame(x = c(mean(contours_plot[contours_plot$level == 250 & contours_plot$tot.tests == 5000,"Ds"]),
                          mean(contours_plot[contours_plot$level == 1000 & contours_plot$tot.tests == 5000,"Ds"])),
                    y = c(5075,5075),
                    label = c("250 infection threshold", "1000 infection threshold"))
  
  segment = data.frame(x = c(max(contours_plot %>% filter(level == 250, d == 0, sens == 0.9) %>% pull(Ds)),
                             max(contours_plot %>% filter(level == 1000, d == 0, sens == 0.9) %>% pull(Ds))),
                       y = 5075/10000, 
                       xend = c(max(contours_plot %>% filter(level == 250, d == 0.3, sens == 1) %>% pull(Ds)),
                                max(contours_plot %>% filter(level == 1000, d == 0.2, sens == 1) %>% pull(Ds))),
                       yend = 5075/10000)
  
  p2 = ggplot()+
    geom_path(data = contours_plot %>% filter(level == 250),
              aes(x = Ds, y = tot.tests/10000, color = d, linetype = as.factor(sens), group = paste(d, as.factor(sens))), 
              size = 1.2)+
    geom_path(data = contours_plot %>% filter(level == 1000),
              aes(x = Ds, y = tot.tests/10000, color = d, linetype = as.factor(sens), group = paste(d, as.factor(sens))), 
              size = 1.2)+
    geom_segment(data = segment, aes(x = x, y = y, xend = xend, yend = yend))+
    geom_label(data = text, aes(x = x, y = y/10000, label = label), fill = "white", label.size = 0)+
    guides(size = FALSE)+
    labs(shape = "Tolerance")+
    scale_color_distiller(palette = "Reds", direction = 0, name = "NPI Intensity", labels = percent_format(accuracy = 1), limits = c(0,0.55))+
    scale_linetype_manual(values = c("dashed", "solid"), labels = c("90%","100%"), name = "Effective\nSensitivity")+
    #scale_size_manual(values = c(1.55, 0.9, 0.6))+ß
    scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)", limits = c(0,8))+
    scale_y_continuous(expand = c(0,0), name = "Test administration\n (% of population tested daily)", labels = percent_format(accuracy = 1), limits = c(0,0.5125))+
    theme_classic(base_size = 13)+
    theme(legend.key.width=unit(1.2,"cm"))
  
  ggsave(file.path(output_path),p2,width = 10, height = 7, units = "in")
}
