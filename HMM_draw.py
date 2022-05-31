class HMM_draw:

    def __init__(self):
        self.radius=1
        self.scale=0.25
        self.text_size = 20*2*self.radius*self.scale
        self.x_origin = 2
        self.y_origin = 2
        self.vertical_spacing = 4
        self.horizontal_spacing = 4
        self.number_of_nodes = 10
        self.figure_size_x = 16
        self.figure_size_y = 9
        self.fig = plt.figure(figsize=(16, 9))
        self.ax = self.fig.add_axes((0, 0, 1, 1))
        self.ax.set_xlim(0, 16)
        self.ax.set_ylim(0, 9)
        #ax.set_facecolor(bg_color)

        self.ax.tick_params(bottom=False, top=False,
                       left=False, right=False)
        self.ax.tick_params(labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)
        return

    def draw_circle(self,x,y):
        theta = np.linspace(0, 2 * np.pi, 100)
        print(np.sin(theta))
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

    def draw_text(self,x,y,text, h_align, v_align, font_size, col):
        text_kwargs = dict(ha=h_align, va=v_align, fontsize=font_size, fontname='Arial', color=col)
        plt.text(x, y, text, **text_kwargs)
        return

    def draw_arrow(self,x1,y1,x2,y2):
        self.ax.annotate("",(x1,y1),(x2,y2),arrowprops=dict(arrowstyle="<|-"))
        return

    def draw_main_block_diamond_arrows(self,x,y):
        self.draw_insert_block_arrows(self.scale*(x*self.horizontal_spacing-self.radius)+self.x_origin, self.scale*(y*self.horizontal_spacing)+self.y_origin)
        #draw to next square
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.y_origin + self.scale * (0 * self.vertical_spacing+self.radius))
        # draw to next delete circle
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing + 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.scale*(2*self.vertical_spacing-0.5*self.radius)+self.y_origin)
        return

    def draw_last_block_diamond_arrows(self, x, y):
        self.draw_insert_block_arrows(self.scale * (x * self.horizontal_spacing - self.radius) + self.x_origin,
                                      self.scale * (y * self.horizontal_spacing) + self.y_origin)
        # draw to next square
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.5 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.5 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing - self.radius), self.y_origin + self.scale * (0 * self.vertical_spacing+self.radius))
        return

    def draw_main_block_delete_arrows(self,x,y):
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing+ 0.707 * self.radius), self.y_origin+self.scale*(y*self.vertical_spacing-0.707*self.radius),
                        self.x_origin + self.scale * ((x+1) * self.horizontal_spacing-0.5*self.radius), self.scale * (0 * self.vertical_spacing+self.radius)+self.y_origin)
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing), self.y_origin+self.scale*(y*self.vertical_spacing-self.radius),
                        self.x_origin + self.scale*x*self.horizontal_spacing, self.y_origin+self.scale*((y-1)*self.vertical_spacing+self.radius))
        self.draw_arrow(self.x_origin + self.scale*(x*self.horizontal_spacing+self.radius), self.y_origin+self.scale*(y*self.vertical_spacing),
                        self.x_origin + self.scale*((x+1)*self.horizontal_spacing-self.radius), self.y_origin+self.scale*(y*self.vertical_spacing))
        return

    def draw_end_block_delete_arrows(self, x, y):
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing + 0.707 * self.radius),
                        self.y_origin + self.scale * (y * self.vertical_spacing - 0.707 * self.radius),
                        self.x_origin + self.scale * ((x + 1) * self.horizontal_spacing-0.5*self.radius), self.scale * (0 * self.vertical_spacing+self.radius)+self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * (x * self.horizontal_spacing),
                        self.y_origin + self.scale * (y * self.vertical_spacing - self.radius),
                        self.x_origin + self.scale * x * self.horizontal_spacing,
                        self.y_origin + self.scale * ((y - 1) * self.vertical_spacing + self.radius))
        return

    def draw_main_match_block_arrows(self,x,y):
        self.draw_arrow(self.x_origin + self.scale * (x+self.radius), self.y_origin,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - self.radius), self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * x, self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * x,
                        self.y_origin + self.scale * (self.vertical_spacing - self.radius))
        self.draw_arrow(self.x_origin+self.scale*(x+self.radius), self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - 0.707 * self.radius),
                        self.y_origin + self.scale * (2 * self.vertical_spacing - 0.707 * self.radius))
        return

    def draw_last_main_match_block_arrows(self, x, y):
        self.draw_arrow(self.x_origin + self.scale * (x+self.radius), self.y_origin,
                        self.x_origin + self.scale * (x+self.horizontal_spacing - self.radius), self.y_origin)
        self.draw_arrow(self.x_origin + self.scale * x, self.y_origin + self.scale * self.radius,
                        self.x_origin + self.scale * x,
                        self.y_origin + self.scale * (self.vertical_spacing - self.radius))
        return

    def draw_insert_block_arrows(self,x,y):
        self.draw_circular_arrow(x, y)
        return

    def draw_circular_arrow(self,x,y):
        #fig = plt.figure(figsize=(9, 9))
        #ax = plt.gca()
        x_radius = self.scale
        y_radius=self.scale
        base_angle = 0 # the 'reference' placement from where angles are measured is 0
        start_angle = 55 # the angle will start at 225 degrees, from a 45 degree angle
        end_angle =247 # the length of the arrow is 270 degrees from the start_angle
        arrow_capstyle = 'round'
        arrow_linestyle = '-'
        arrow_linewidth = 1
        arrow_colour = 'black'
        # ========Line
        arc = Arc([x, y], x_radius, y_radius, angle=start_angle, theta1=base_angle, theta2=end_angle, capstyle=arrow_capstyle, linestyle=arrow_linestyle, lw=arrow_linewidth, color=arrow_colour)
        self.ax.add_patch(arc)
        # ========Create the arrow head
        arrow_head_x = x + (x_radius / 2) * np.cos(rad(end_angle + start_angle))  # Do trig to determine end position
        arrow_head_y = y + (y_radius / 2) * np.sin(rad(end_angle + start_angle))
        # Create triangle as arrow head # (x,y) # number of vertices # radius # orientation
        self.ax.add_patch(RegularPolygon((arrow_head_x, arrow_head_y), 3, x_radius/15, rad(start_angle+end_angle), color=arrow_colour))
        return

    def draw(self):
        self.draw_begin()
        for i in range(self.number_of_nodes):
            self.draw_square((i+1), 0)
            self.draw_text(self.x_origin + self.scale * ((i+1) * self.horizontal_spacing),
                               self.y_origin + self.scale * (0 * self.vertical_spacing), 'M' + str(i+1), 'center',
                               'center', self.text_size,
                               'black')
            if i < self.number_of_nodes-1:
                self.draw_main_match_block_arrows((i+1)*self.horizontal_spacing, 0)
            else:
                self.draw_last_main_match_block_arrows((i+1)*self.horizontal_spacing, 0)
            self.draw_diamond((i + 1), 1)
            if i < self.number_of_nodes-1:
                self.draw_main_block_diamond_arrows((i+1), 1)
            else:
                self.draw_last_block_diamond_arrows(i+1, 1)
            self.draw_circle((i + 1), 2)
            self.draw_text(self.x_origin + self.scale * ((i + 1) * self.horizontal_spacing),
                           self.y_origin + self.scale * (2 * self.vertical_spacing), 'D' + str(i + 1), 'center',
                           'center', self.text_size,
                           'black')
            if i < self.number_of_nodes - 1:
                self.draw_main_block_delete_arrows(i+1, 2)
            else:
                self.draw_end_block_delete_arrows(i+1, 2)
            # if i < self.number_of_nodes-1:
            #      self.draw_main_block_delete_arrows(i, 2)
            # else:
            #     self.draw_end_block_delete_arrows(i, 2)
        self.draw_end()
        return

    def draw_figure_title(self, text):
        self.draw_text((self.scale*(self.number_of_nodes+2)*self.horizontal_spacing+self.x_origin)*0.5, self.scale*3 * self.vertical_spacing+self.y_origin, 'My HMM 1', 'center', 'center', self.text_size, 'black')
        return