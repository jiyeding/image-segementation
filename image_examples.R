library(imager)
library(RSpectra)

### elephant example ####
elephantcimg <- load.image("elephant.jpg")
elephantcimg <- imresize(elephantcimg, 1/10)
elephantcimg <- grayscale(elephantcimg)

plot(elephantcimg)

elephant <- as.matrix(elephantcimg)
data <- elephant

# no landmarks#
system.time(image.njw(data=data,b=25,sigma_p = 0.02,x=3))
#landmarks
system.time(image.lsc(data=data,B=25,b=3,sigma_p = 0.03,x=3))

##### tree example ####
treescimg <- load.image("trees.jpg")
treescimg <- imresize(treescimg, 1/50)
treescimg <- grayscale(treescimg)

plot(treescimg)

trees <- as.matrix(treescimg)
data <- trees

# no landmarks #
system.time(image.njw(data=data,b=11,sigma_p = 0.05,x=4))
# landmarks #
system.time(image.lsc(data=data,B=25,b=3,sigma_p = 0.03,x=4))

##### moutain example ####
mountaincimg <- load.image("mountain.jpg")
mountaincimg <- imresize(mountaincimg, 1/30)
mountaincimg <- grayscale(mountaincimg)

plot(mountaincimg)

mountain <- as.matrix(mountaincimg)
data <- mountain

# no landmark #
system.time(image.njw(data=data,b=11,sigma_p = 0.1,x=3))
# landmarks #
system.time(image.lsc(data=data,B=17,b=3,sigma_p = 0.1,x=3)) 