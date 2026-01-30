import pygame
import numpy as np
import math

windows_size = (800, 800)
class Interface:
    def polygon_points(self ,center , size):
        height = size
        width = math.sqrt(3) * size
        n0 = (center[0] - height, center[1] + width)
        n1 = (center[0] + height, center[1] + width)
        n2 = (center[0] + 2 * height, center[1])
        n3 = (center[0] + height, center[1] - width)
        n4 = (center[0] - height, center[1] - width)
        n5 = (center[0] - 2 * height, center[1])
        return [n0, n1, n2, n3, n4, n5]


    def draw_polygon(self,center,color):
        pygame.draw.polygon(
            self.screen,
            color,
            self.polygon_points(center, self.cell_size)
        )
        pygame.draw.polygon(
            self.screen,
            (255, 255, 255),
            self.polygon_points(center, self.cell_size),
            width=1
        )

    def draw_grid(self):
        self.screen.fill((20, 20, 20))
        width = math.sqrt(3) * self.cell_size
        height = self.cell_size
        initial_center = [2*self.cell_size,2*self.cell_size]
        center = initial_center.copy()
        for i in range(self.matrix.shape[1]):
            j=0
            self.draw_polygon(center, self.color[self.matrix[i, j]])
            for j in range(self.matrix.shape[0]-1):
                if j % 2 == 0:
                    center[0] += 3 * height
                    center[1] += width
                    self.draw_polygon([center[0], center[1]], self.color[self.matrix[i, j]])
                else:
                    center[0] += 3 * height
                    center[1] -= width
                    self.draw_polygon([center[0], center[1]], self.color[self.matrix[i, j]])
            center[0] = initial_center[0]
            center[1] = i * 2 * width + initial_center[1]
        text = self.font.render(f'Iteration {self.iteration_counter}', True, (255,255,255))
        textRect = text.get_rect()
        textRect.bottomright = self.windows_size
        self.screen.blit(text, textRect)
        pygame.display.flip()






    def __init__(self,windows_size, controler):
        self.windows_size = windows_size
        self.controler=controler
        pygame.init()
        self.screen = pygame.display.set_mode(windows_size)
        pygame.display.set_caption("Lizard")
        running = False
        self.matrix = self.controler.getGrid()
        self.color = [(0,0,0), (255, 255, 255)]
        self.cell_size = ((windows_size[0]) / (3.5 * self.matrix.shape[0]))*0.98
        clock = pygame.time.Clock()
        self.iteration_counter = 0
        self.matrix_history = []
        self.font = pygame.font.Font('freesansbold.ttf', 24)

        self.matrix_history.append(self.matrix)

        # Main loop
        while True:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: #close
                    pygame.quit()
                    return

                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_SPACE: # Start/ stop
                        running = not running
                    elif event.key == pygame.K_RIGHT:
                        self.matrix = step(self.matrix)
                        self.iteration_counter += 1
                        self.matrix_history.append(self.matrix)
                    elif event.key == pygame.K_LEFT:
                        self.iteration_counter -= 1
                        self.matrix = self.matrix_history[self.iteration_counter]
                        self.matrix_history = self.matrix_history[:-1]
                    elif event.key == pygame.K_c: # Clear the board
                        self.matrix.fill(0)

            
            if running:
                self.matrix = self.controler.getGrid()
                self.iteration_counter += 1
                self.matrix_history.append(self.matrix)
            
            self.draw_grid()
            clock.tick(60)
