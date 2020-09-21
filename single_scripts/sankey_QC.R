# Library
library(networkD3)
library(dplyr)
library(htmlwidgets)

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
    source=c("group_A","group_A", "group_B", "group_C", "group_C", "group_E"), 
    target=c("group_C","group_D", "group_E", "group_F", "group_G", "group_H"), 
    value=c(2,3, 2, 3, 1, 3)
)
links <- data.frame(
    source = c("All",  "All",   "Dead", "Dead", "PEEK",    "PEEK",    "PEEK",     "PEEK",          "Drift A",   "Drift B",   "No drift",  "Corrected",    "Corrected",  "Night method", "Night method", "Day method", "Day method"), 
    target = c("Dead", "Alive", "PEEK", "Al",   "Drift A", "Drift B", "No drift", "Not corrected", "Corrected", "Corrected", "Corrected", "Night method", "Day method", "QC 1",         "QC 2",         "QC 1",       "QC 2"), 
    value1 =  c(26,     7,       21,     5,      11,        2,         3,          5,               11,          2,           3,            11,            5,            10,             1,              1,            4),
    value2 =  c(11,     3,       10,     1,      5,         0,         2,          3,               5,           0,           2,            6,             1,            2,              4,              1,            0),
    value3 =  c(22,     7,       17,     5,      6,         4,         2,          5,               6,           4,           2,            9,             3,            6,              3,              2,            1),
    value4 =  c(10,     2,       10,     0,      4,         1,         2,          3,               4,           1,           2,            2,             5,            2,              0,              4,            1)
)

links$value = links$value1 + links$value2 + links$value3 + links$value4

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
    name = c(as.character(links$source), 
           as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

levels(nodes$name)[1] = paste0("Aluminium (", links$value[4], ")")
levels(nodes$name)[2] = paste0("Alive (", links$value[2], ")")
levels(nodes$name)[3] = paste0("All (", links$value[1] + links$value[2], ")")
levels(nodes$name)[4] = paste0("Corrected (", links$value[12] + links$value[13], ")")
levels(nodes$name)[5] = paste0("Day method (", links$value[13], ")")
levels(nodes$name)[6] = paste0("Dead (", links$value[1], ")")
levels(nodes$name)[7] = paste0("Drift A (", links$value[5], ")")
levels(nodes$name)[8] = paste0("Drift B (", links$value[6], ")")
levels(nodes$name)[9] = paste0("Night method (", links$value[12], ")")
levels(nodes$name)[10] = paste0("No drift (", links$value[7], ")")
levels(nodes$name)[11] = paste0("Not corrected (", links$value[8], ")")
levels(nodes$name)[12] = paste0("PEEK (", links$value[3], ")")
levels(nodes$name)[13] = paste0("QC 1 (", links$value[14] + links$value[16], ")")
levels(nodes$name)[14] = paste0("QC 2 (", links$value[15] + links$value[17], ")")

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, fontSize=13)
p
saveWidget(p, file="sankey_QC.html", selfcontained=FALSE)
