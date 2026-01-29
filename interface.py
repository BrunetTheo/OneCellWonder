import pygame
import numpy as np
import math
from dataclasses import dataclass

@dataclass
class Config: # Basic param for the screen / game
    cell_size = 10
    grid_w = 30
    grid_h = 30
    debug = False

def polygon_points(center, size):
    height = size
    width = math.sqrt(3) * size
    n0 = (center[0] - height, center[1] + width)
    n1 = (center[0] + height, center[1] + width)
    n2 = (center[0] + 2 * height, center[1])
    n3 = (center[0] + height, center[1] - width)
    n4 = (center[0] - height, center[1] - width)
    n5 = (center[0] - 2 * height, center[1])
    return [n0, n1, n2, n3, n4, n5]


def main():
    cfg = Config()
    pygame.init()
    screen = pygame.display.set_mode((cfg.grid_w * 3*cfg.cell_size + 2*cfg.cell_size, cfg.grid_h * 2 * math.sqrt(3)*cfg.cell_size))
    pygame.display.set_caption("my beautiful window that will hopefully work one day")
    running = False
    matrix = np.random.randint(0,2,size=(cfg.grid_w, cfg.grid_h))
    color = [(0,0,0), (255, 255, 255)]

    def debug(txt):
        if cfg.debug:
            print(txt)

    def draw_polygon(center, color):
        pygame.draw.polygon(
            screen,
            color,
            polygon_points(center, cfg.cell_size)
        )
        pygame.draw.polygon(
            screen,
            (255, 255, 255),
            polygon_points(center, cfg.cell_size),
            width=1
        )

    def draw_grid():
        screen.fill((20, 20, 20))
        width = math.sqrt(3) * cfg.cell_size
        height = cfg.cell_size
        initial_center = [25,25]
        center = initial_center.copy()
        for i in range(matrix.shape[1]):
            j=0
            draw_polygon(center, color[matrix[i, j]])
            for j in range(matrix.shape[0]-1):
                if j % 2 == 0:
                    center[0] += 3 * height
                    center[1] += width
                    draw_polygon([center[0], center[1]], color[matrix[i, j]])
                else:
                    center[0] += 3 * height
                    center[1] -= width
                    draw_polygon([center[0], center[1]], color[matrix[i, j]])
            center[0] = initial_center[0]
            center[1] = i * 2 * width + initial_center[1]
        pygame.display.flip()

    # Main loop
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: #close
                pygame.quit()
                return
        draw_grid()


if __name__ == "__main__":
    main()
