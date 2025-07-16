# This file indicates that colors and shapes of metadata across various
# let's keep the visualizations consistent! :) 


# water color scheme for taxa
ch4_colors <- c("Methanobacteriales" = "#2E624D",
                "Methanocellales" = "#558D7B", 
                "Methanofastidiosales" =  "#148F77",
                "Methanomassiliicoccales" = "#009E73",
                "Methanomicrobiales" = "#43BA8F",
                "Methanosarcinales_A_2632" =  "#48C9B0",
                "Methanomethylicales" = "#90CEC0",
                "Methanotrichales" = "#A3E4D7",
                "Methylococcales"= "#7D3560",
                "Methylomirabilales" = "#CC79A7",
                "Methylacidiphilales"= "#EFB6D6"
)



solar_colors <- c(
  "Open" = "#67a9cf",  #"#76A7CB",
  "In Progress" = "#A6B7C6",
  "FPV" = "#ef8a62" # "#C07A5B"
)

solar_colorsv1.2 <- c(
  "No Solar" = "#76A7CB",
  "In Progress" = "#A6B7C6",
  "Solar" = "#C07A5B"
)
solar_colorsv1 <- c(
  "No Solar" = "#B3C493",
  "In Progress" = "#A6B7C6",
  "Solar" = "#005373"
)

solar_colors1 <- c(
  "No Solar" = "#BABB74", 
  "In Progress" = "#2A9D8F", 
  "Solar" = "#264653"
)

solar_colors2 <- c(
  "Solar" = "#440154FF",
  "In Progress" = "#21908CFF", 
  "No Solar" = "#FDE725FF"
)

solar_colors3 <- c(
  "Solar" = "darkolivegreen3", 
  "In Progress" = "darkgoldenrod1", 
  "No Solar" = "cadetblue2"
)


solar_shapes <- c(
  "Solar" = 12, 
  "In Progress" = 3,
  "No Solar" = 17)

month_shapes <- c(
  "August" = 23, 
  "September" = 21,
  "October" = 22,
  "Control" = 24
)

fraction_colors <- c(
  "Particle" = "firebrick3", 
  "Free" = "goldenrod1", 
  "Whole" = "darkorange2",
  "Control" = "grey")

month_colors <- c(
  "October" = "dodgerblue4",
  "September" = "dodgerblue2",
  "August" = "#00ADA7",
  "Control" = "grey")

pond_colors <- c(
  "123" = "#440154FF",
  "124" = "#DB8300",
  "125" = "#6A2626",
  "128" = "#F2C621",
  "131" = "#9FDA3AFF",
  "132" = "#E3A5D6")


pond_colors1 <- c(
  "133" = "navy",
  "131" = "darkslategray2",
  "127" = "gold",
  "123" = "plum2",
  "134" = "red1",
  "132" = "darkorange",
  "128" = "firebrick",
  "124" = "limegreen",
  "135" = "magenta",
  "125" = "royalblue")

pond_colors2 <- c(
  "123" = "#440154FF",
  "124" = "#414487FF",
  "125" = "#2A788EFF",
  "128" = "#22A884FF",
  "131" = "#7AD151FF",
  "132" = "#FDE725FF"
)

pond_colors3 <- c(
  "123" = "#4AC16DFF",
  "124" = "#FDE725FF",
  "125" = "#9FDA3AFF",
  "128" = "#365C8DFF",
  "131" = "#46337EFF",
  "132" = "#440154FF"
)

pond_shapes <- c(
  "131" = 15,
  "123" = 10,
  "132" = 17,
  "128" = 16, #18,
  "124" = 12,
  "125" = 13)

depth_colors <- c(
  "B" = "#3A4332",
  "S" = "#D48792",
  "SD" = "#774C1E")

sample_colors <- c(
  "Water" = "#7BAEA0",
  "Sediment" = "#774C1E"
)

sample_shapes <- c(
  "Water" = 17,
  "Sediment" = 22
)


depth_shapes <- c(
  "B" = 19, 
  "S" = 2,
  "SD" = 22
)

year_shapes <- c(
  "2022" = 14, 
  "2023" = 17,
  "2024"= 19
)

year_colors <- c(
  "2022" = "#FDE725FF", 
  "2023" = "#821554",
  "2024" = "#D48349"
)


# Set the ggplot theme
theme_set(theme_bw() + 
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 9)))

date_collected_colors <- c(
  "2023-06-12" = "#4B0049",
  "2023-06-29" = "#5D014F",
  "2023-07-20" = "#700853",
  "2023-08-11" = "#821554",
  "2023-08-22" = "#932252",
  "2023-09-13" = "#963D4E", 
  "2023-10-03" = "#985350",
  "2023-10-25" = "#95675B", 
  "2024-06-20" = "#B16D51",
  "2024-07-11" = "#E48F43",
  "2024-08-21" = "#F29C3B",
  "2024-09-11" = "#FFBA54"
)

# Set the phylum colors
phylum_colors <- c(
  Acidobacteriota = "navy", 
  Actinobacteriota = "darkslategray2", 
  Armatimonadota = "deeppink1",
  Alphaproteobacteria = "plum2", 
  Bacteroidota = "gold", 
  Betaproteobacteria = "plum1", 
  Bdellovibrionota = "red1",
  Chloroflexi="black", 
  Crenarchaeota = "firebrick",
  Cyanobacteria = "limegreen",
  Deltaproteobacteria = "grey", 
  Desulfobacterota="magenta",
  Euryarchaeota = "blue",
  Firmicutes = "#3E9B96",
  Gammaproteobacteria = "greenyellow",
  "Marinimicrobia (SAR406 clade)" = "yellow",
  Myxococcota = "#B5D6AA",
  Nitrospirota = "palevioletred1",
  Patescibacteria = "darkmagenta",
  Proteobacteria = "royalblue",
  Planctomycetota = "darkorange", 
  "SAR324 clade(Marine group B)" = "olivedrab",
  Thermoplasmatota = "green",
  Verrucomicrobiota = "darkorchid1", 
  Campylobacterota = "darkorange", 
  Spirochaetota = "forestgreen", 
  Nanoarchaeota = "lightpink1",
  Dependentiae = "greenyellow", 
  MBNT15 = "darkorange3", 
  Micrarchaeota = "darkseagreen3", 
  Sumerlaeota = "deeppink4", 
  Sva0485 = "goldenrod3", 
  "WPS-2" = "darkolivegreen2", 
  Gemmatimonadota = "coral2",
  Other = "darkgrey",
  Euryarchaeota = "#6089B5",
  Halobacterota = "#CD622E",
  Nanoarchaeota = "#D76E9A",
  Sva0485 = "thistle")

# Other = "grey")
# 
# # 
# Acadia = list(c("#212E52", "#444E7E", "#8087AA", "#B7ABBC", "#F9ECE8", "#FCC893", "#FEB424", "#FD8700", "#D8511D"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=TRUE),
# Arches = list(c("#1A3D82", "#0C62AF", "#4499F5", "#8FCAFD", "#F2F2F2", "#F0AC7D", "#CD622E", "#B14311", "#832B0F"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=TRUE),
# Arches2 = list(c("#3A1F46", "#7F4B89", "#B46DB3", "#E3A5D6", "#F3DAE4"), c(1, 2, 3, 4, 5), colorblind=TRUE),
# Banff = list(c("#006475", "#00A1B7", "#55CFD8", "#586028", "#898928", "#616571", "#9DA7BF"), c(2, 5, 1, 6, 3, 7, 4), colorblind=FALSE),
# BryceCanyon = list(c("#882314", "#C0532B", "#CF932C", "#674D53", "#8C86A0", "#724438", "#D5AB85"), c(1, 5, 2, 7, 4, 3, 6), colorblind=FALSE),
# CapitolReef = list(c("#291919", "#532A34", "#7C5467", "#878195", "#AEB2B7", "#D4D9DD"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
# Chamonix = list(c("#008FF8", "#B6AA0D", "#E2C2A2", "#E23B0E", "#F2C621", "#196689"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# CraterLake = list(c("#1D4A79", "#794C23", "#6B7444", "#6089B5", "#BF9785", "#275E4D", "#807B7F"), c(1, 2, 3, 4, 5, 6, 7), colorblind=FALSE),
# Cuyahoga = list(c("#E07529", "#FAAE32", "#7F7991", "#A84A00", "#5D4F36", "#B39085"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
# DeathValley = list(c("#8C2B0E", "#C5692D", "#FEB359", "#132F5B", "#435F90", "#68434E", "#B47E83"), c(1, 5, 7, 2, 6, 3, 4), colorblind=TRUE),
# Denali = list(c("#20223E", "#3F3F7B", "#278192", "#00B089", "#2EEA8C", "#8FF7BD"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# Everglades = list(c("#345023", "#596C0B", "#83A102", "#003B68", "#426F86", "#7A712F"), c(3, 4, 1, 6, 5, 2), colorblind=FALSE),
# Glacier = list(c("#01353D", "#088096", "#58B3C7", "#7AD4E4", "#B8FCFC"), c(1, 2, 3, 4, 5), colorblind=TRUE),
# GrandCanyon = list(c("#521E0F", "#9C593E", "#DDA569", "#3F4330", "#8E7E3C", "#2A4866", "#6592B0"), c(2, 6, 3, 4, 7, 1, 5), colorblind=FALSE),
# Halekala = list(c("#722710", "#A3844D", "#675243", "#A85017", "#838BAA"), c(1, 2, 3, 4, 5), colorblind=TRUE),
# IguazuFalls = list(c("#415521", "#97AD3D", "#4C3425", "#7F6552", "#5A8093", "#9FBAD3"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# KingsCanyon = list(c("#613921", "#A77652", "#F2C27B", "#AAC9ED", "#44637D", "#8E949F"), c(1, 5, 6, 3, 2, 4), colorblind=TRUE),
# LakeNakuru = list(c("#D76E9A", "#A1ACC8", "#AD3C36", "#332627", "#EACACF", "#AA6B77"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# Olympic = list(c("#3A4330", "#426737", "#75871B", "#BAB97D", "#FAF3CE", "#FDE16A", "#F9B40E", "#E88C23", "#A25933"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=FALSE),
# Redwood = list(c("#5E3B49", "#9B5F6B", "#BA817D", "#325731", "#6A9741", "#5F4E2F"), c(2, 5, 6, 3, 4, 1), colorblind=FALSE),
# RockyMtn = list(c("#274C31", "#A3AEB5", "#2F4B6A", "#8F8081", "#3F7156", "#6F89A7", "#5B5443"), c(1, 2, 3, 4, 5, 6, 7), colorblind=FALSE),
# Saguaro = list(c("#127088", "#C85729", "#92874B", "#CD8A39", "#AC3414", "#57643C"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# SmokyMtns = list(c("#42511A", "#889D35", "#D3D175", "#B50200", "#DA6C41", "#7C6E66", "#BCAFA6"), c(1, 4, 2, 6, 3, 5, 7), colorblind=FALSE),
# SouthDowns = list(c("#948D2A", "#D5B44D", "#89A4BF", "#F1D6B6", "#9B8358", "#577291"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# Torres = list(c("#2F397A", "#7391BD", "#894846", "#E9988C", "#535260", "#B7A7A6", "#785838", "#C68D61", "#4F6008", "#93995C"), c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), colorblind=FALSE),
# Triglav = list(c("#386EC2", "#B5B5B2", "#990006", "#625D0A", "#B9741F", "#213958"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
# WindCave = list(c("#2F100E", "#6C3322", "#B07159", "#C9A197", "#E0CDCD"), c(1, 2, 3, 4, 5), colorblind=TRUE),
# Volcanoes = list(c("#082544", "#1E547D", "#79668C", "#DE3C37", "#F2DC7E"), c(1, 2, 3, 4, 5), colorblind=TRUE),
# # Yellowstone = list(c("#0067A2", "#DFCB91", "#CB7223", "#289A84", "#7FA4C2", "#AF7E56"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
# Yosemite = list(c("#293633", "#3D5051", "#6B7F7F", "#87A1C7", "#516B95", "#304F7D"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE)
# 
# rainier = c(lake="#465177", ragwort="#E4C22B", lodge="#965127",
#             trees="#29483A", ground="#759C44", winter_sky="#9FB6DA",
#             paintbrush="#DF3383"),
# washington_pass = c(trees="#31543B", stone="#48628A", tips="#94AA3D",
#                     road="#7F9CB1", sunbreak="#D9D1BE", stump="#3E3C3A"),
# palouse = c(snake="#2D3F4A", wheat="#C0A43D", fallow="#8A6172",
#             hills="#748A52", canyon="#CCBA98", sky="#69A2E4"),
# forest = c(trees="#254029", stream="#1E212F", fern="#516F25",
#            bark="#3A270A", mountains="#40666F"),
# larch = c(larch="#D2A554", shrub="#626B5D", rock="#8C8F9E",
#           moss="#858753", sky="#A4BADF", dirt="#D3BEAF"),
# coast = c(surf="#7BAEA0", sea="#386276", rocks="#3A4332", sand="#7A7D6F",
#           sunset="#D9B96E", sky="#BED4F0"),
# san_juan = c(trees="#21281D", grass="#CA884C", sea="#3A5775",
#              driftwood="#BAAF9F", clouds="#C9DCE2"),
# uw = c(purple="#483778", gold="#D2C28B", brick="#7E4837", cherry="#D48792",
#        stone="#6B7471"),
# fort_worden = c(sea="#4D5370", shrub="#263C19", lighthouse="#A6B7C6",
#                 rust="#8F6D3F"),
# skagit = c(red="#B51829", yellow="#EFC519", violet="#831285",
#            orange="#CC7810", purple="#C886E8", mountains="#3E6E94"),
# # flag = c(green="#247F5B", yellow="#E2CD70", blue="#5DB3DB", tan="#FFDCCD"),
# # # CONTINUOUS
# sound_sunset = c("#001E36", "#0F2649", "#352C5A", "#533369", "#6E3D71",
#                  "#814C74", "#925C78", "#A26D7C", "#BC7A7D", "#D08B79",
#                  "#DE9F71", "#E7B565", "#EBCC5C", "#E7E55C", "#DCFF6C"),
# ferries = c("#241E33", "#1E2C54", "#003F69", "#005373", "#006372",
#             "#007273", "#008174", "#349075", "#579E78", "#77AB7E",
#             "#97B788", "#B3C493", "#CFD28F", "#F2DE83", "#FFE66C"),
# forest_fire = c("#2C0915", "#40111D", "#551A23", "#6A2626", "#7E3225",
#                 "#92401D", "#A64F00", "#B96000", "#CB7100", "#DB8300",
#                 "#E69800", "#EDAD00", "#F3C307", "#FAD753", "#FFEC7A"),
# sea = c("#1C1B23", "#222541", "#1F325A", "#094270", "#005282", "#0B6184",
#         "#2F6F89", "#457D90", "#598C99", "#6B9AA3", "#7CA9AE", "#8DB9B9",
#         "#9DC8C5", "#ADD8D2", "#BEE8DF"),
# sea_star = c("#4B0049", "#5D014F", "#700853", "#821554", "#932252",
#              "#963D4E", "#985350", "#95675B", "#B16D51", "#C3774D",
#              "#D48349", "#E48F43", "#F29C3B", "#FFAA3B", "#FFBA54"),
# volcano = c("#29272C", "#323335", "#3E3F3F", "#4A4B4B", "#555857",
#             "#626561", "#71726B", "#837E77", "#8C8D87", "#989B98",
#             "#A7A8A8", "#B5B7B7", "#C2C5C8", "#D3D3D9", "#E7E0E8"),
# baker = c("#627F9A", "#7684A7", "#8B88B1", "#A08CB9", "#B490BF", "#C795C3",
#           "#D69CC3", "#E1A5C0", "#E9AFBE", "#EEBABD", "#F2C6BE", "#F4D2C2",
#           "#F4DEC9", "#F5E9D3", "#F6F4E5"),
# diablo = c("#172512", "#0E320C", "#003F0F", "#004C1A", "#005929",
#            "#00663C", "#007051", "#007B63", "#008575", "#279085",
#            "#3E9B96", "#51A6A6", "#63B2B7", "#75BDC7", "#87C8D8"),
# puget = c("#1D3024", "#123C2E", "#00483E", "#005352", "#005C67", "#01657D",
#           "#386B91", "#5F6FA3", "#8272B1", "#A175B8", "#B87DB8", "#CB86B6",
#           "#DA92B3", "#E69FAF", "#EFAEAB"),
# mountains = c("#002733", "#003141", "#003B50", "#004661", "#005172",
#               "#005C84", "#006797", "#0071AA", "#3A7BB5", "#5E85BE",
#               "#7990C7", "#8F9CCF", "#A4A8D8", "#B7B4E1", "#C8C2EA"),
# gorge = c("#00352A", "#00402F", "#004C33", "#005736", "#076338", "#236E38",
#           "#387937", "#4D8434", "#658E3C", "#7C984D", "#91A15E", "#A5AB70",
#           "#B7B583", "#C8C095", "#D9CBA7"),
# foothills = c("#001C28", "#00262D", "#003233", "#003E39", "#004B3E",
#               "#005841", "#006542", "#00723F", "#007F38", "#18893C",
#               "#529249", "#739C5B", "#8FA572", "#A6AF8C", "#B9BAAF"),
# footbridge = c("#32221C", "#42291C", "#50311D", "#5E391D", "#6B431C",
#                "#774C1E", "#7D5933", "#846545", "#8B7155", "#927C64",
#                "#9A8873", "#A39481", "#ACA08E", "#B6AC9B", "#C1B8A7"),
# olympic = c("#2A2A40", "#3E3E63", "#535285", "#6B6B98", "#8484A4",
#             "#9E9EB4", "#B9B9C5", "#D4D4DA", "#B3BDB0", "#93A58E",
#             "#738F6A", "#517942", "#396128", "#2D4823", "#232F1F"),
# lopez = c("#A3730F", "#B2833A", "#C19456", "#D0A570", "#DEB689", "#EBC8A3",
#           "#F6DBBF", "#F1F1F1", "#DBDCFF", "#C9CAFA", "#B7B8F0", "#A5A7E5",
#           "#9496DA", "#8385CE", "#7275C2"),
# vantage = c("#473065", "#5B497F", "#70629B", "#857EAB", "#949EB5",
#             "#A4BDC6", "#C0DAD6", "#EFF2E5", "#DCD797", "#C5BB61",
#             "#AD9F3C", "#9A8224", "#89641E", "#774718", "#632B14"),
# stuart = c("#4F3173", "#6F457F", "#8D5C8C", "#A87499", "#C18FA8",
#            "#D8AAB8", "#EDC7CB", "#FFE2DE", "#F6C7AE", "#D8B07A",
#            "#AF9C57", "#85884B", "#617344", "#425C3D", "#2C4534")
# )
