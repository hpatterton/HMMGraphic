import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad

class HMM_draw:

    def __init__(self, number_of_nodes=5, transition_values_array = 0):
        self.radius=1
        self.scale=0.3
        self.text_size = 20*2*self.radius*self.scale
        self.x_origin = 2
        self.y_origin = 2
        self.vertical_spacing = 4
        self.horizontal_spacing = 4
        self.number_of_nodes = number_of_nodes
        self.figure_size_x = 16
        self.figure_size_y = 9
        self.fig = plt.figure(figsize=(16, 9))
        self.ax = self.fig.add_axes((0, 0, 1, 1))
        self.ax.set_xlim(0, 16)
        self.ax.set_ylim(0, 9)
        self.arrow_text_coordinates = {'MD': (.8, 1.30), 'MI': (0.12, 0.45), 'MM': (0.46, 0.1), 'ID': (0.45, 1.7),
                                       'II': (-0.27, 0.75),
                                       'IM': (0.46, 0.35), 'DD': (0.5, 2.1), 'DI': (0.12, 1.5), 'DM': (0.33, 1.1)}
        self.number_of_conserved_columns = 10
        self.dictionary_of_transition_paths = {'MM': 0, 'MI': 1, 'MD': 2, 'IM': 3, 'II': 4, 'ID': 5, 'DM': 6, 'DI': 7, 'DD': 8}
        if transition_values_array == 0:
            self.transition_value_array = [[0.123] * len(self.dictionary_of_transition_paths) for i in range(self.number_of_conserved_columns + 1)]
        else:
            self.transition_value_array = transition_values_array
        #ax.set_facecolor(bg_color)


        self.ax.tick_params(bottom=False, top=False,
                       left=False, right=False)
        self.ax.tick_params(labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
        return

    def draw_circle(self,x,y):
        theta = np.linspace(0, 2 * np.pi, 100)
        self.ax.plot(self.scale*(x*self.horizontal_spacing + self.radius * np.cos(theta))+self.x_origin, self.scale*(y*self.vertical_spacing + self.radius * np.sin(theta))+self.y_origin,color="midnightblue",)
        return

    def draw_square(self,x,y):
        self.ax.plot((self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin,self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin),(self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin,self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin),color="midnightblue",linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin), color="midnightblue" ,lw=1,linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing + self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin), color="midnightblue" ,linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing + self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing - self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing - self.radius)+self.y_origin), color="midnightblue" ,linestyle='-')
        return

    def draw_diamond(self,x,y):
        self.ax.plot((self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing)+self.x_origin), (self.scale*(y*self.vertical_spacing)+self.y_origin, self.scale*(y*self.vertical_spacing+self.radius)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing)+self.x_origin, self.scale*(x*self.horizontal_spacing+self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing+self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing)+self.y_origin), color="midnightblue", lw=1,
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing+self.radius)+self.x_origin, self.scale*(x*self.horizontal_spacing)+self.x_origin), (self.scale*(y*self.vertical_spacing)+self.y_origin, self.scale*(y*self.vertical_spacing-self.radius)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.ax.plot((self.scale*(x*self.horizontal_spacing)+self.x_origin, self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin), (self.scale*(y*self.vertical_spacing-self.radius)+self.y_origin, self.scale*(y*self.vertical_spacing)+self.y_origin), color="midnightblue",
                     linestyle='-')
        self.draw_text(self.x_origin + self.scale * (x * self.horizontal_spacing),
                       self.y_origin + self.scale * (y * self.vertical_spacing), 'I' + str(x), 'center', 'center', self.text_size,
                       'black')

        return

    def show_drawing(self):
        plt.show()
        return

    def draw_begin(self):
        # draw 'Begin'
        x=0
        y=0
        self.draw_square(x, y)
        self.draw_main_match_block_arrows(x,y)
        self.draw_text(self.x_origin, self.y_origin, 'Begin', 'center', 'center', self.text_size, 'black')

        # draw insert diamond
        self.draw_diamond(x, y+1)
        self.draw_main_block_diamond_arrows(x, y+1)
        return

    def draw_end(self):
        self.draw_square((self.number_of_nodes+1), 0)
        self.draw_text(self.x_origin + self.scale * ((self.number_of_nodes+1) * self.horizontal_spacing),
                       self.y_origin + self.scale * (0 * self.vertical_spacing), 'End', 'center',
                       'center', self.text_size,
                       'black')
        return

    def draw_transition_values(self, transition_value_array, dictionary_of_transition_paths = {'MM': 0, 'MI': 1, 'MD': 2, 'IM': 3, 'II': 4, 'ID': 5, 'DM': 6, 'DI': 7, 'DD': 8}):
        number_of_nodes = len(transition_value_array)
        list_of_keys = list(dictionary_of_transition_paths.keys())
        number_of_keys = len(list_of_keys)
        for node in range(number_of_nodes):
            for key in list_of_keys:
                self.write_text_at_transition_arrow(node, key, transition_value_array[node][dictionary_of_transition_paths[key]], 'black', 10)
        return


    def draw_text(self,x,y,text, h_align, v_align, font_size, col):
        text_kwargs = dict(ha=h_align, va=v_align, fontsize=font_size, fontname='Arial', color=col)
        plt.text(x, y, text, **text_kwargs)
        return

    def write_text_at_transition_arrow(self, node, transition_key, text, text_color='black', text_size=10):
        number = float(text)
        if number == 0:
            text = '0'
        elif number == 1:
            text = '1.0'
        else:
            text = ('%5.3f' % float(text)).rstrip('0').rstrip('.')
        start_key = ['MD', 'MI', 'MM', 'ID', 'II', 'IM']
        body_key = ['MD', 'MI', 'MM', 'ID', 'II', 'IM', 'DD', 'DI', 'DM']
        end_key = ['MI', 'MM', 'II', 'IM', 'DI', 'DM']
        if node == 0 and transition_key in start_key:
            self.draw_text(self.x_origin + self.scale * (
                        (node + self.arrow_text_coordinates[transition_key][0]) * self.horizontal_spacing),
                           self.y_origin + self.scale * (
                                       self.arrow_text_coordinates[transition_key][1] * self.vertical_spacing), text,
                           'center', 'center', text_size, text_color)
        elif node > 0 and node < self.number_of_nodes and transition_key in body_key:
            self.draw_text(self.x_origin + self.scale * (
                        (node + self.arrow_text_coordinates[transition_key][0]) * self.horizontal_spacing),
                           self.y_origin + self.scale * (
                                       self.arrow_text_coordinates[transition_key][1] * self.vertical_spacing), text,
                           'center', 'center', text_size, text_color)

        elif node == self.number_of_nodes and transition_key in end_key:
            self.draw_text(self.x_origin + self.scale * (
                        (node + self.arrow_text_coordinates[transition_key][0]) * self.horizontal_spacing),
                           self.y_origin + self.scale * (
                                       self.arrow_text_coordinates[transition_key][1] * self.vertical_spacing), text,
                           'center', 'center', text_size, text_color)
        return

    def draw_transmission_text(self):

        # start
        keys = ['MD', 'MI', 'MM', 'ID', 'II', 'IM']
        node = 0
        for i in range(len(keys)):
            self.draw_text(self.x_origin + self.scale * ((node + self.arrow_text_coordinates[keys[i]][0]) * self.horizontal_spacing),
                           self.y_origin + self.scale * (self.arrow_text_coordinates[keys[i]][1] * self.vertical_spacing), str(keys[i]), 'center',
                           'center', self.text_size,
                           'blue')
        # main body
        keys = ['MD', 'MI', 'MM', 'ID', 'II', 'IM', 'DD', 'DI', 'DM']
        for node in range(1, self.number_of_nodes):
            for i in range(len(keys)):
                self.draw_text(self.x_origin + self.scale * ((node + self.arrow_text_coordinates[keys[i]][0]) * self.horizontal_spacing),
                       self.y_origin + self.scale * (self.arrow_text_coordinates[keys[i]][1] * self.vertical_spacing), str(keys[i]), 'center',
                       'center', self.text_size,
                       'blue')
        # end
        keys = ['MI', 'MM', 'II', 'IM', 'DI', 'DM']
        node = self.number_of_nodes
        for i in range(len(keys)):
            self.draw_text(self.x_origin + self.scale * ((node + self.arrow_text_coordinates[keys[i]][0]) * self.horizontal_spacing),
                           self.y_origin + self.scale * (self.arrow_text_coordinates[keys[i]][1] * self.vertical_spacing), str(keys[i]),
                           'center',
                           'center', self.text_size,
                           'blue')

        return

    def draw_arrow(self,x1,y1,x2,y2,arrow_color='black',arrow_width=1):
        self.ax.annotate("",(x1,y1),(x2,y2),arrowprops=dict(arrowstyle="<|-", color=arrow_color, lw=arrow_width))
        return

    def draw_main_match_block_arrows(self,x,y):
        self.draw_arrow_special(x, y, 'MM', 'black', 1)
        self.draw_arrow_special(x, y, 'MI', 'black', 1)
        self.draw_arrow_special(x, y, 'MD', 'black', 1)
        return

    def draw_last_main_match_block_arrows(self, x, y):
        self.draw_arrow_special(x, y, 'MM', 'black', 1)
        self.draw_arrow_special(x, y, 'MI', 'black', 1)
        return

    def draw_main_block_diamond_arrows(self,x,y):
        self.draw_arrow_special(x, y, 'II', 'black', 1)
        self.draw_arrow_special(x, y, 'IM', 'black', 1)
        self.draw_arrow_special(x, y, 'ID', 'black', 1)
        return

    def draw_last_block_diamond_arrows(self, x, y):
        self.draw_arrow_special(x, y, 'II', 'black', 1)
        self.draw_arrow_special(x, y, 'IM', 'black', 1)

        return

    def draw_main_block_delete_arrows(self,x,y):
        self.draw_arrow_special(x, y, 'DM', 'black', 1)
        self.draw_arrow_special(x, y, 'DI', 'black', 1)
        self.draw_arrow_special(x, y, 'DD', 'black', 1)
        return

    def draw_end_block_delete_arrows(self, x, y):
        self.draw_arrow_special(x, y, 'DM', 'black', 1)
        self.draw_arrow_special(x, y, 'DI', 'black', 1)
        return

    def draw_insert_block_arrows(self,x,y, arrow_color, arrow_width):
        self.draw_circular_arrow(x, y, arrow_color, arrow_width)
        return

    def draw_arrow_special(self, x, y, key, arrow_color='black', arrow_width=1):
        if key == 'DM':
            self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.707 * self.radius),
                            self.y_origin + self.scale * (y * self.vertical_spacing - 0.707 * self.radius),
                            self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - 0.5 * self.radius),
                            self.scale * (0 * self.vertical_spacing + self.radius) + self.y_origin, arrow_color, arrow_width)
        if key == 'DI':
            self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing),
                            self.y_origin + self.scale * (y * self.vertical_spacing - self.radius),
                            self.x_origin + self.scale * x * self.horizontal_spacing,
                            self.y_origin + self.scale * ((y - 1) * self.vertical_spacing + self.radius), arrow_color, arrow_width)
        if key == 'DD':
            self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + self.radius),
                            self.y_origin + self.scale * (y * self.vertical_spacing),
                            self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius),
                            self.y_origin + self.scale * (y * self.vertical_spacing), arrow_color, arrow_width)
        if key == 'MM':
            self.draw_arrow(self.x_origin + self.scale * (x + self.radius), self.y_origin,
                            self.x_origin + self.scale * (x + self.horizontal_spacing - self.radius), self.y_origin, arrow_color, arrow_width)
        if key == 'MI':
            self.draw_arrow(self.x_origin + self.scale * x, self.y_origin + self.scale * self.radius,
                            self.x_origin + self.scale * x,
                            self.y_origin + self.scale * (self.vertical_spacing - self.radius), arrow_color, arrow_width)
        if key == 'MD':
            self.draw_arrow(self.x_origin + self.scale * (x + self.radius), self.y_origin + self.scale * self.radius,
                            self.x_origin + self.scale * (x + self.horizontal_spacing - 0.707 * self.radius),
                            self.y_origin + self.scale * (2 * self.vertical_spacing - 0.707 * self.radius), arrow_color, arrow_width)
        if key == 'II':
            self.draw_insert_block_arrows(self.scale * (x * self.horizontal_spacing - self.radius) + self.x_origin,
                                          self.scale * (y * self.horizontal_spacing) + self.y_origin, arrow_color, arrow_width)
        if key == 'IM':
            self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                            self.y_origin + self.scale * (y * self.vertical_spacing - 0.5 * self.radius),
                            self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius),
                            self.y_origin + self.scale * (0 * self.vertical_spacing + self.radius), arrow_color, arrow_width)
        if key == 'ID':
            self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                            self.y_origin + self.scale * (y * self.vertical_spacing + 0.5 * self.radius),
                            self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius),
                            self.scale * (2 * self.vertical_spacing - 0.5 * self.radius) + self.y_origin, arrow_color, arrow_width)

        return

    def draw_circular_arrow(self, x, y, arrow_color, arrow_width):
        x_radius = self.scale
        y_radius=self.scale
        base_angle = 0 # the 'reference' placement from where angles are measured is 0
        start_angle = 55 # the angle will start at 225 degrees, from a 45 degree angle
        end_angle =247 # the length of the arrow is 270 degrees from the start_angle
        arrow_capstyle = 'round'
        arrow_linestyle = '-'
        # ========Line
        arc = Arc([x, y], x_radius, y_radius, angle=start_angle, theta1=base_angle, theta2=end_angle, capstyle=arrow_capstyle, linestyle=arrow_linestyle, lw=arrow_width, color=arrow_color)
        self.ax.add_patch(arc)
        # ========Create the arrow head
        arrow_head_x = x + (x_radius / 2) * np.cos(rad(end_angle + start_angle))  # Do trig to determine end position
        arrow_head_y = y + (y_radius / 2) * np.sin(rad(end_angle + start_angle))
        # Create triangle as arrow head # (x,y) # number of vertices # radius # orientation
        self.ax.add_patch(RegularPolygon((arrow_head_x, arrow_head_y), 3, x_radius/15, rad(start_angle+end_angle), color=arrow_color, lw=arrow_width))
        return

    def draw_figure_title(self, text):
        self.draw_text((self.scale*(self.number_of_nodes+1)*self.horizontal_spacing+self.x_origin)*0.5, self.scale*3 * self.vertical_spacing+self.y_origin, 'My HMM 1', 'center', 'center', self.text_size, 'black')
        return

    def draw(self):
        self.draw_begin()
        for i in range(1, self.number_of_nodes + 1):
            self.draw_square(i, 0)
            self.draw_text(self.x_origin + self.scale * (i * self.horizontal_spacing),
                           self.y_origin + self.scale * (0 * self.vertical_spacing), 'M' + str(i), 'center',
                           'center', self.text_size,
                           'black')
            if i < self.number_of_nodes:
                self.draw_main_match_block_arrows(i * self.horizontal_spacing, 0)
            else:
                self.draw_last_main_match_block_arrows(i * self.horizontal_spacing, 0)
            self.draw_diamond(i, 1)
            if i < self.number_of_nodes:
                self.draw_main_block_diamond_arrows(i, 1)
            else:
                self.draw_last_block_diamond_arrows(i, 1)
            self.draw_circle(i, 2)
            self.draw_text(self.x_origin + self.scale * (i * self.horizontal_spacing),
                           self.y_origin + self.scale * (2 * self.vertical_spacing), 'D' + str(i), 'center',
                           'center', self.text_size,
                           'black')
            if i < self.number_of_nodes:
                self.draw_main_block_delete_arrows(i, 2)
            else:
                self.draw_end_block_delete_arrows(i, 2)
        self.draw_transition_values(self.transition_value_array)
        #self.draw_transmission_text()
        self.draw_end()

        return