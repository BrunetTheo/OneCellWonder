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

    def point_in_polygon(self, point, polygon):
        """Check if a point is inside a polygon using ray casting algorithm"""
        x, y = point
        n = len(polygon)
        inside = False
        p1x, p1y = polygon[0]
        for i in range(n + 1):
            p2x, p2y = polygon[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def get_cell_at_mouse(self, mouse_pos):
        """Find which cell the mouse is hovering over"""
        width = math.sqrt(3) * self.cell_size
        height = self.cell_size
        initial_center = [2*self.cell_size, 2*self.cell_size]
        center = initial_center.copy()
        
        for i in range(self.matrix.shape[1]):
            j = 0
            polygon = self.polygon_points(center, self.cell_size)
            if self.point_in_polygon(mouse_pos, polygon):
                return (i, j)
            
            for j in range(self.matrix.shape[0]-1):
                if j % 2 == 0:
                    center[0] += 3 * height
                    center[1] += width
                else:
                    center[0] += 3 * height
                    center[1] -= width
                
                polygon = self.polygon_points(center, self.cell_size)
                if self.point_in_polygon(mouse_pos, polygon):
                    return (i, j)
            
            center[0] = initial_center[0]
            center[1] = i * 2 * width + initial_center[1]
        
        return None

    def draw_tooltip(self, cell_coords):
        """Draw a tooltip showing the first 10 gene contents"""
        if cell_coords is None:
            return
        
        i, j = cell_coords
        gene_content = self.controler.cellGrid.gene_content[i, j]
        
        # Get first 10 genes
        genes_to_show = min(10, len(gene_content))
        
        # Create tooltip background
        tooltip_width = 130
        tooltip_height = 30 + genes_to_show * 25
        mouse_pos = pygame.mouse.get_pos()
        
        # Position tooltip near mouse, but keep it on screen
        tooltip_x = mouse_pos[0] + 15
        tooltip_y = mouse_pos[1] + 15
        
        if tooltip_x + tooltip_width > self.windows_size[0]:
            tooltip_x = mouse_pos[0] - tooltip_width - 15
        if tooltip_y + tooltip_height > self.windows_size[1]:
            tooltip_y = mouse_pos[1] - tooltip_height - 15
        
        # Draw semi-transparent background
        tooltip_surface = pygame.Surface((tooltip_width, tooltip_height))
        tooltip_surface.set_alpha(230)
        tooltip_surface.fill((40, 40, 40))
        pygame.draw.rect(tooltip_surface, (255, 255, 255), 
                        tooltip_surface.get_rect(), 2)
        self.screen.blit(tooltip_surface, (tooltip_x, tooltip_y))
        
        # Draw title
        title_text = self.small_font.render(f'Cell ({i}, {j}) ', 
                                           True, (255, 255, 100))
        self.screen.blit(title_text, (tooltip_x + 10, tooltip_y + 5))
        
        # Draw gene values
        for idx in range(genes_to_show):
            gene_text = self.small_font.render(
                f'Gene {idx}: {int(gene_content[idx])}', 
                True, (255, 255, 255))
            self.screen.blit(gene_text, 
                           (tooltip_x + 10, tooltip_y + 30 + idx * 25))

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
        
        # Draw iteration counter
        text = self.font.render(f'Iteration {self.iteration_counter}', True, (255,255,255))
        textRect = text.get_rect()
        textRect.bottomright = self.windows_size
        self.screen.blit(text, textRect)
        
        # Draw tooltip if hovering over a cell
        mouse_pos = pygame.mouse.get_pos()
        hovered_cell = self.get_cell_at_mouse(mouse_pos)
        if hovered_cell is not None:
            self.draw_tooltip(hovered_cell)
        
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
        self.small_font = pygame.font.Font('freesansbold.ttf', 18)  # For tooltip

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
                        self.controler.update()
                        self.matrix = self.controler.getGrid()
                        self.iteration_counter += 1
                        self.matrix_history.append(self.matrix)
                    elif event.key == pygame.K_LEFT:
                        self.iteration_counter -= 1
                        self.matrix = self.matrix_history[self.iteration_counter]
                        self.matrix_history = self.matrix_history[:-1]
                    elif event.key == pygame.K_a: # Clear the board
                        self.controler.show = 0
                        self.matrix = self.controler.getGrid()
                    elif event.key == pygame.K_z: # Clear the board
                        self.controler.show = 1
                        self.matrix = self.controler.getGrid()
                    elif event.key == pygame.K_r: # Clear the board
                        self.controler.show = -1
                        self.matrix = self.controler.getGrid()
            
            if running:
                self.controler.update()
                self.matrix = self.controler.getGrid()
                self.iteration_counter += 1
                self.matrix_history.append(self.matrix)
            
            self.draw_grid()
            clock.tick(60)