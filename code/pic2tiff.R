library(magick)

img <- image_read("PNAS/figs/Figure 3.png")
img <- image_convert(img, format = "tiff")

image_write(
  img,
  path = "Figure 3.tiff",
  format = "tiff",
  density = 600,      
  compression = "lzw" 
)
