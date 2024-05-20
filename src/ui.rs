use macroquad::miniquad::window::screen_size;
use macroquad::prelude::*;
use crate::FIELD_SIZE;

use crate::state::State;

pub struct UI {
    pub day_color: Color,
    pub night_color: Color,
    pub background_color: Color
}

impl UI {
    pub fn draw(&self, state: &State) {
        let (w, h) = screen_size();
        let size = std::cmp::min_by(w, h, |a, b| a.partial_cmp(&b).unwrap()) * 0.6;
        let size = size - size % FIELD_SIZE as f32;
        let cell_size = size / FIELD_SIZE as f32;

        let left = (w - size) / 2.;
        let top = (h - size) / 2.;
        let font_size = size * 0.05;

        clear_background(self.background_color);

        self.draw_field(state, left, top, cell_size);
        self.draw_stats(state, left, top + size * 1.01, size, font_size);
    }

    fn draw_field(&self, state: &State, left: f32, top: f32, cell_size: f32) {
        for y in 0..state.field.grid.len() {
            let row = &state.field.grid[y];
            for x in 0..row.len() {
                let color = if row[x] { self.day_color } else { self.night_color };
                draw_rectangle(left + x as f32 * cell_size, top + y as f32 * cell_size, cell_size, cell_size, color)
            }
        }
        let ball_radius = cell_size / 2.;
        draw_circle(left + state.day_ball.pos.x * cell_size, top + state.day_ball.pos.y * cell_size, ball_radius, self.night_color);
        draw_circle(left + state.night_ball.pos.x * cell_size, top + state.night_ball.pos.y * cell_size, ball_radius, self.day_color);
    }

    fn draw_stats(&self, state: &State, left: f32, top: f32, width: f32, font_size: f32) {
        let text = &format!("day {} | night {}", state.field.day_count, state.field.night_count);
        let text_dimensions = measure_text(text, None, font_size as u16, 1.);
        let x = left + (width - text_dimensions.width) / 2.;
        let y = top;
        draw_text(text, x, y + text_dimensions.offset_y, font_size, BLACK);
    }
}
