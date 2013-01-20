
science <- read.csv("sf.csv", header=FALSE)
colnames(science) <- c('name',
                       'WATCH BIRDS',
                       'READ BOOKS ON ANIMALS',
                       'READ BOOKS ON PLANTS',
                       'WATCH GRASS CHANGE',
                       'FIND BOTTLES AND CANS',
                       'LOOK UP STRANGE ANIMAL OR PLANT',
                       'WATCH ANIMAL MOVE',
                       'LOOK IN SIDEWALK CRACKS',
                       'LEARN WEED NAMES',
                       'LISTEN TO BIRD SING',
                       'FIND WHERE ANIMAL LIVES',
                       'GO TO MUSEUM',
                       'GROW GARDEN',
                       'LOOK AT PICTURES OF PLANTS',
                       'READ ANIMAL STORIES',
                       'MAKE A MAP',
                       'WATCH WHAT ANIMALS EAT',
                       'GO ON PICNIC',
                       'GO TO ZOO',
                       'WATCH BUGS',
                       'WATCH BIRD MAKE NEST',
                       'FIND OUT WHAT ANIMALS EAT',
                       'WATCH A RAT',
                       'FIND OUT WHAT FLOWERS LIVE ON',
                       'TALK W FRIENDS ABOUT PLANTS')
dim(science)
for (cx in 2:dim(science)[2]) {
  science[[cx]] <- ordered(science[[cx]], levels=0:2, labels=c('dislike','neutral','like'))
}

colnames(science.people) <- c('trait','se')

science.people <- as.data.frame(science.people)
for (c in colnames(science)) {
  science.people[[c]] <- science[[c]]
}

save(science.people, file="science.people.rda")
save(science.items, file="science.items.rda")

sfif <- read.csv("SFIF.txt")
sfpf <- read.csv("SFPF.txt")
sfxf <- read.csv("SFXF.txt")
