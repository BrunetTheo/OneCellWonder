import pygame
import numpy as np
import math

windows_size = (800, 800)

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
    pygame.init()
    screen = pygame.display.set_mode(windows_size)
    pygame.display.set_caption("my beautiful window that will hopefully work one day")
    running = False
    matrix = np.random.randint(0,2,size=(500, 500))
    color = [(0,0,0), (255, 255, 255)]
    cell_size = (windows_size[0] - 50) / (3.5 * matrix.shape[0])
    clock = pygame.time.Clock()
    iteration_counter = 0
    matrix_history = []
    font = pygame.font.Font('freesansbold.ttf', 24)

    def draw_polygon(center, color):
        pygame.draw.polygon(
            screen,
            color,
            polygon_points(center, cell_size)
        )
        pygame.draw.polygon(
            screen,
            (255, 255, 255),
            polygon_points(center, cell_size),
            width=1
        )

    def draw_grid():
        screen.fill((20, 20, 20))
        width = math.sqrt(3) * cell_size
        height = cell_size
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
        text = font.render(f'Iteration {iteration_counter}', True, (255,255,255))
        textRect = text.get_rect()
        textRect.bottomright = windows_size
        screen.blit(text, textRect)
        pygame.display.flip()
    matrix_history.append(matrix)

    # Main loop
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT: #close
                pygame.quit()
                return

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE: # Start/ stop
                    running = not running
                if event.key == pygame.K_RIGHT:
                    matrix = step(matrix)
                    iteration_counter += 1
                    matrix_history.append(matrix)
                if event.key == pygame.K_LEFT:
                    iteration_counter -= 1
                    matrix = matrix_history[iteration_counter]
                    matrix_history = matrix_history[:-1]

        
        if running:
            matrix = step(matrix)
            iteration_counter += 1
            matrix_history.append(matrix)
        
        draw_grid()

if __name__ == "__main__":
    main()
