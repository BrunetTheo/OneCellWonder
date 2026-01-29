import pygame
import numpy as np
import math
from dataclasses import dataclass

@dataclass
class Config: # Basic param for the screen / game
    cell_size = 10
    grid_w = 100
    grid_h = 100
    debug = False

def polygon_points(center, size):
    height = size
    width = math.sqrt(3) * size
    n0 = (center, center + 2 * height)
    n1 = (center + width, center + height)
    n2 = (center + width, center - height)
    n3 = (center, center - 2 * height)
    n4 = (center - width, center - height)
    n5 = (center - width, center + height)
    return [n0, n1, n2, n3, n4, n5]


def main():
    cfg = Config()
    pygame.init()
    screen = pygame.display.set_mode((cfg.grid_w * cfg.cell_size, cfg.grid_h * cfg.cell_size))
    pygame.display.set_caption("my beautiful window that will hopefully work one day")
    running = False

    def debug(txt):
        if cfg.debug:
            print(txt)

    def draw():
        screen.fill((20, 20, 20))
        pygame.draw.polygon(screen, (255, 255, 255), polygon_points(500, 50))
        pygame.display.flip()
    # Main loop
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: #close
                pygame.quit()
                return
        draw()


if __name__ == "__main__":
    main()
