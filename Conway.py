import pygame
import numpy as np
from dataclasses import dataclass


def step(grid: np.ndarray) -> np.ndarray:
    neighbors = sum(
        np.roll(np.roll(grid, dy, axis=0), dx, axis=1)
        for dy in (-1, 0, 1)
        for dx in (-1, 0, 1)
        if not (dx == 0 and dy == 0)
    )

    birth = (neighbors == 3) & (grid == 0) # Rule 1
    survive = ((neighbors == 2) | (neighbors == 3)) & (grid == 1) # Rule 2
    
    return (birth | survive).astype(np.uint8)

@dataclass
class Config: # Basic param for the screen / game
    cell_size = 10
    grid_w = 100
    grid_h = 100
    fps = 15  # simulation speed when running
    emerg = 0.20 # Proba of the random uni
    debug = False


@dataclass
class Gun: # Class to store the coordinate for the gosper glider gun
    position_b1 = [[4,0],
                   [5,0],
                   [4,1],
                   [5,1]]
    position_b2 = [[2,34],
                   [3,34],
                   [2,35],
                   [3,35]]
    position_o1 = [[4,10],[5,10],[6,10],
                   [3,11],[7,11],
                   [2,12],[8,12],
                   [2,13],[8,13],
                   [5,14],
                   [3,15],[7,15],
                   [4,16],[5,16],[6,16],
                   [5,17]]
    position_o2 = [[2,20],[3,20],[4,20],
                   [2,21],[3,21],[4,21],
                   [1,22],[5,22],
                   [0,24],[1,24],[5,24],[6,24]]
    gun = position_b1+position_b2+position_o1+position_o2




def main():
    cfg = Config()
    gun = Gun()
    pygame.init()
    screen = pygame.display.set_mode((cfg.grid_w * cfg.cell_size, cfg.grid_h * cfg.cell_size))
    pygame.display.set_caption("my beautiful window that will hopefully work one day")
    clock = pygame.time.Clock()
    grid = np.zeros((cfg.grid_h, cfg.grid_w), dtype=np.uint8)
    running = False

    def debug(txt):
        if cfg.debug:
            print(txt)

    def draw():
        screen.fill((20, 20, 20))
        # Draw alive cells
        ys, xs = np.where(grid == 1)
        for y, x in zip(ys, xs):
            rect = pygame.Rect(
                x * cfg.cell_size, y * cfg.cell_size, cfg.cell_size, cfg.cell_size
            )
            pygame.draw.rect(screen, (220, 220, 220), rect)

        for x in range(cfg.grid_w):
            pygame.draw.line(
                screen, (40, 40, 40),
                (x * cfg.cell_size, 0),
                (x * cfg.cell_size, cfg.grid_h * cfg.cell_size),
                1
            )
        for y in range(cfg.grid_h):
            pygame.draw.line(
                screen, (40, 40, 40),
                (0, y * cfg.cell_size),
                (cfg.grid_w * cfg.cell_size, y * cfg.cell_size),
                1
            )
        pygame.display.flip()
    # Main loop
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: #close
                pygame.quit()
                return

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE: # Start/ stop
                    running = not running
                elif event.key == pygame.K_UP: # speed UP
                    cfg.fps = min(120, cfg.fps + 5) # Cap at 120
                elif event.key == pygame.K_DOWN: # speed DOWN
                    cfg.fps = max(1, cfg.fps - 5) # minimum of 1 gen/s

                elif event.key == pygame.K_c: # Clear the board
                    grid.fill(0)
                
                elif event.key == pygame.K_r: # Randomly fill the 15% of the grid
                    grid = (np.random.rand(cfg.grid_h, cfg.grid_w) < cfg.emerg).astype(np.uint8)
                elif event.key == pygame.K_s: # create a slider
                    mx, my = pygame.mouse.get_pos()
                    gx = mx // cfg.cell_size
                    gy = my // cfg.cell_size
                    debug(f"position X = {gx}\nPosition Y:{gy}")
                    if 0 <= gx-1 and gx+1 < cfg.grid_w and 0 <= gy-1 and gy+1 < cfg.grid_h: 
                        grid[gy+1, gx+1] = 1  
                        grid[gy, gx+1] = 1  
                        grid[gy-1, gx+1] = 1  
                        grid[gy-1, gx] = 1  
                        grid[gy, gx-1] = 1    
                elif event.key == pygame.K_g: # create a Gosper glider gun
                    mx, my = pygame.mouse.get_pos()
                    gx = mx // cfg.cell_size
                    gy = my // cfg.cell_size
                    debug(f"position X = {gx}\nPosition Y:{gy}")
                    if 0 <= gx and gx+36 < cfg.grid_w and 0 <= gy and gy+9 < cfg.grid_h: 
                        for pos in gun.gun:
                            grid[gy+ pos[0],gx + pos[1]] = 1
                        


            if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1: # Place one cell
                mx, my = pygame.mouse.get_pos()
                gx = mx // cfg.cell_size
                gy = my // cfg.cell_size
                debug(f"position X = {gx}\nPosition Y:{gy}")
                if 0 <= gx < cfg.grid_w and 0 <= gy < cfg.grid_h:
                    grid[gy, gx] ^= 1

        if running: #update the grid before display
            grid = step(grid)
        draw()
        clock.tick(cfg.fps if running else 60) # set the speed


if __name__ == "__main__":
    main()
