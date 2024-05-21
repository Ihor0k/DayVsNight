use macroquad::prelude::*;

use crate::state::State;
use crate::ui::UI;

mod ui;
mod state;

const FIELD_SIZE: usize = 20;
const CELL_SIZE_PX: u32 = 20;
const DAY_COLOR_HEX: u32 = 0xCBD9D9;
const NIGHT_COLOR_HEX: u32 = 0x16313E;

const BACKGROUND_COLOR_HEX: u32 = 0xE7EBEC;
const SPEED: f32 = 40.;

fn window_conf() -> Conf {
    Conf {
        window_title: "Day vs Night".to_owned(),
        window_height: 800,
        window_width: 800,
        window_resizable: false,
        ..Default::default()
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let time_ms = (miniquad::date::now() * 1000.) as u64;
    rand::srand(time_ms);

    let mut state = State::default();
    let day_color = Color::from_hex(DAY_COLOR_HEX);
    let night_color = Color::from_hex(NIGHT_COLOR_HEX);
    let background_color = Color::from_hex(BACKGROUND_COLOR_HEX);

    let ui = UI {
        day_color,
        night_color,
        background_color,
    };

    loop {
        let mut dt = get_frame_time();

        if dt > 1. {
            dt = dt % 1.
        }
        state.next_step(dt * SPEED);

        ui.draw(&state);
        next_frame().await
    }
}
